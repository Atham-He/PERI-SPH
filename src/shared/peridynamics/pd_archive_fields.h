/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_ARCHIVE_FIELDS_H
#define PD_ARCHIVE_FIELDS_H
#include "pd_canonical_archive.h"
#include <memory>
#include <type_traits>
#include <utility>

namespace SPH::peridynamics::canonical_archive
{
namespace field_detail
{
template<class T,class=void>struct Matrix:std::false_type{};
template<class T>struct Matrix<T,std::void_t<typename T::Scalar,decltype(T::RowsAtCompileTime),decltype(T::ColsAtCompileTime)>>:std::true_type{};

template<class T>std::vector<uint64_t> elementShape()
{
    if constexpr(Matrix<T>::value)
    {
        static_assert(std::is_same_v<typename T::Scalar,double>,"Archive matrix scalars require binary64");
        static_assert(T::RowsAtCompileTime>0 && T::ColsAtCompileTime>0,"Only explicit fixed-size matrix elements are admitted");
        if constexpr(T::ColsAtCompileTime==1)return {uint64_t(T::RowsAtCompileTime)};
        else return {uint64_t(T::RowsAtCompileTime),uint64_t(T::ColsAtCompileTime)};
    }
    else
    {
        static_assert(std::is_integral_v<T> || std::is_same_v<T,double>,"Unsupported archive scalar");
        return {};
    }
}
template<class T>std::vector<uint64_t> arrayShape(uint64_t n)
{
    auto shape=elementShape<T>();shape.insert(shape.begin(),n);elements(shape);return shape;
}
template<class T>uint64_t integer(T value)
{
    static_assert(std::is_integral_v<T> && sizeof(T)<=8,"Integer profile");
    if constexpr(std::is_signed_v<T>)
        if(value<0)throw std::invalid_argument("archive fields: negative unsigned-profile integer");
    return uint64_t(value);
}
template<class T>T fromInteger(uint64_t value)
{
    static_assert(std::is_integral_v<T> && sizeof(T)<=8,"Integer profile");
    if(value>uint64_t(std::numeric_limits<T>::max()))throw std::invalid_argument("archive fields: integer narrowing");
    return T(value);
}
template<class T>void appendReals(std::vector<double> &values,const T &value)
{
    if constexpr(Matrix<T>::value)
        for(int i=0;i<T::RowsAtCompileTime;++i)for(int j=0;j<T::ColsAtCompileTime;++j)values.push_back(value(i,j));
    else
    {
        static_assert(std::is_same_v<T,double>,"Real profile");values.push_back(value);
    }
}
template<class T>T takeReals(const std::vector<double> &values,size_t &offset)
{
    if constexpr(Matrix<T>::value)
    {
        T value;
        for(int i=0;i<T::RowsAtCompileTime;++i)for(int j=0;j<T::ColsAtCompileTime;++j)value(i,j)=values.at(offset++);
        return value;
    }
    else
    {
        static_assert(std::is_same_v<T,double>,"Real profile");return values.at(offset++);
    }
}
}

/** Typed registration adapter. Explicit row-major components, never Eigen or
 * vector object storage. Named scopes keep the control and dense map distinct.
 * This layer preserves fields, but cannot qualify model semantics or restart.
 */
class FieldWriter
{
    Builder &archive_;std::string prefix_;
    std::string key(const std::string &name)const{return prefix_+name;}
  public:
    explicit FieldWriter(Builder &archive,std::string prefix={}):archive_(archive),prefix_(std::move(prefix)){}
    FieldWriter scope(const std::string &name)const{return FieldWriter(archive_,key(name)+"/");}
    template<class T>void addScalar(const std::string &name,const T &value)
    {
        const auto shape=field_detail::elementShape<T>();
        if constexpr(std::is_same_v<T,uint8_t>)archive_.bytes(key(name),shape,{value});
        else if constexpr(std::is_integral_v<T>)archive_.integers(key(name),shape,{field_detail::integer(value)});
        else
        {
            std::vector<double> data;field_detail::appendReals(data,value);archive_.reals(key(name),shape,data);
        }
    }
    template<class Array>void addArray(const std::string &name,const Array &values)
    {
        using T=typename Array::value_type;
        const auto shape=field_detail::arrayShape<T>(values.size());
        archive_.integers(key(name)+"/count",{}, {uint64_t(values.size())});
        if constexpr(std::is_same_v<T,uint8_t>)archive_.bytes(key(name),shape,Bytes(values.begin(),values.end()));
        else if constexpr(std::is_integral_v<T>)
        {
            std::vector<uint64_t> data;data.reserve(values.size());
            for(T x:values)data.push_back(field_detail::integer(x));
            archive_.integers(key(name),shape,data);
        }
        else
        {
            std::vector<double> data;data.reserve(elements(shape));
            for(const auto &x:values)field_detail::appendReals(data,x);
            archive_.reals(key(name),shape,data);
        }
    }
};

class FieldReader
{
    const Candidate &archive_;std::string prefix_;
    std::shared_ptr<std::set<std::string>> names_;
    std::string key(const std::string &name)const
    {
        std::string result=prefix_+name;
        if(!names_->insert(result).second)throw std::invalid_argument("archive fields: duplicate schema registration");
        return result;
    }
    FieldReader(const Candidate &archive,std::string prefix,std::shared_ptr<std::set<std::string>> names)
        :archive_(archive),prefix_(std::move(prefix)),names_(std::move(names)){}
  public:
    explicit FieldReader(const Candidate &archive):archive_(archive),names_(std::make_shared<std::set<std::string>>()){}
    FieldReader scope(const std::string &name)const{return FieldReader(archive_,prefix_+name+"/",names_);}
    void finish()const{archive_.requireExactNames(*names_);}
    template<class T>void addScalar(const std::string &name,T &value)
    {
        const auto label=key(name);
        const auto shape=field_detail::elementShape<T>();
        if constexpr(std::is_same_v<T,uint8_t>)value=archive_.bytes(label,shape).at(0);
        else if constexpr(std::is_integral_v<T>)value=field_detail::fromInteger<T>(archive_.integers(label,shape).at(0));
        else
        {
            const auto data=archive_.reals(label,shape);size_t offset=0;value=field_detail::takeReals<T>(data,offset);
        }
    }
    template<class Array>void addArray(const std::string &name,Array &values)
    {
        using T=typename Array::value_type;
        const auto count=archive_.integers(key(name+"/count"),{}).at(0);
        const auto shape=field_detail::arrayShape<T>(count);const auto label=key(name);
        // Validate payload type and dimensions before allocating the receiver.
        if constexpr(std::is_same_v<T,uint8_t>)
        {
            const auto data=archive_.bytes(label,shape);Array temporary(data.begin(),data.end());values=std::move(temporary);
        }
        else if constexpr(std::is_integral_v<T>)
        {
            const auto data=archive_.integers(label,shape);Array temporary;temporary.reserve(size_t(count));
            for(uint64_t x:data)temporary.push_back(field_detail::fromInteger<T>(x));
            values=std::move(temporary);
        }
        else
        {
            const auto data=archive_.reals(label,shape);Array temporary;temporary.reserve(size_t(count));size_t offset=0;
            for(uint64_t i=0;i<count;++i)temporary.push_back(field_detail::takeReals<T>(data,offset));
            values=std::move(temporary);
        }
    }
};

/** There is intentionally no existing/live receiver parameter. Decode and
 * validate a newly owned candidate; callers may commit it only afterward.
 * Component registration is not a complete session or output transaction.
 */
template<class T,class Validate>T decodeDetached(const Candidate &archive,Validate validate)
{
    T temporary;FieldReader reader(archive);temporary.registerCheckpoint(reader);reader.finish();
    validate(std::as_const(temporary));return temporary;
}
}
#endif
