/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_CANONICAL_ARCHIVE_H
#define PD_CANONICAL_ARCHIVE_H

// Standalone byte layer. Not yet bound to a PD session or a live filesystem
// transaction. It never writes through live simulation pointers.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include "pd_sha256.h"

namespace SPH::peridynamics::canonical_archive
{
using Bytes=std::vector<uint8_t>;
static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,
              "Canonical mechanical scalars require IEEE binary64");
inline constexpr uint64_t maximum_bytes=16*1024*1024,maximum_fields=256,
    maximum_elements=1024*1024,maximum_name=128;
inline constexpr std::array<uint8_t,8> magic{'P','E','R','I','A','R','0','1'},
    footer{'P','E','R','I','E','N','D','1'};

/** Use the installed kernel SHA256 implementation, not a hand-rolled hash.
 * No alternative provider or algorithm fallback is permitted in this profile.
 * The provider has a separate known-vector and Python hashlib assessment.
 */
inline Digest sha256(const Bytes &bytes,size_t count)
{
    if(count>bytes.size())throw std::invalid_argument("canonical archive: digest size");
    Sha256Stream stream;stream.update(bytes.data(),count);return stream.finish();
}

inline void appendU64(Bytes &out,uint64_t x)
{for(size_t k=0;k<8;++k)out.push_back(uint8_t(x>>(8*k)));}
inline uint64_t doubleBits(double x)
{uint64_t bits;std::memcpy(&bits,&x,8);return bits;}
inline double fromBits(uint64_t bits)
{double x;std::memcpy(&x,&bits,8);return x;}
inline void validName(const std::string &name)
{
    if(name.empty() || name.size()>maximum_name)
        throw std::invalid_argument("canonical archive: field name size");
    for(unsigned char c:name)
        if(!((c>='A'&&c<='Z')||(c>='a'&&c<='z')||(c>='0'&&c<='9')||c=='_'||c=='/'||c=='.'||c=='-'))
            throw std::invalid_argument("canonical archive: field name alphabet");
}
inline uint64_t elements(const std::vector<uint64_t> &dimensions)
{
    if(dimensions.size()>3)throw std::invalid_argument("canonical archive: rank");
    uint64_t count=1;
    for(uint64_t n:dimensions)
    {
        if(n>maximum_elements)throw std::invalid_argument("canonical archive: dimension bound");
        if(n==0){count=0;continue;}
        if(count>maximum_elements/n)throw std::invalid_argument("canonical archive: element bound");
        count*=n;
    }
    return count;
}
enum class Type:uint8_t{Byte=1,U64=2,F64=3};
inline size_t width(Type type)
{
    if(type==Type::Byte)return 1;
    if(type==Type::U64 || type==Type::F64)return 8;
    throw std::invalid_argument("canonical archive: type");
}
struct Field
{
    Type type;
    std::vector<uint64_t> dimensions;
    Bytes data;
};

class Builder
{
    uint64_t generation_;
    uint64_t payload_bytes_=32;
    std::map<std::string,Field> fields_;
    void add(std::string name,Type type,std::vector<uint64_t> shape,Bytes data)
    {
        validName(name);const uint64_t count=elements(shape);
        if(data.size()!=count*width(type))throw std::invalid_argument("canonical archive: field size");
        const uint64_t increment=18+name.size()+8*shape.size()+data.size();
        if(increment>maximum_bytes-40-payload_bytes_)throw std::invalid_argument("canonical archive: payload byte bound");
        if(fields_.size()>=maximum_fields || !fields_.emplace(std::move(name),Field{type,std::move(shape),std::move(data)}).second)
            throw std::invalid_argument("canonical archive: field count or duplicate");
        payload_bytes_+=increment;
    }
  public:
    explicit Builder(uint64_t explicit_generation):generation_(explicit_generation)
    {if(!generation_)throw std::invalid_argument("canonical archive: explicit generation required");}
    void bytes(std::string name,std::vector<uint64_t> shape,Bytes values)
    {add(std::move(name),Type::Byte,std::move(shape),std::move(values));}
    void integers(std::string name,std::vector<uint64_t> shape,const std::vector<uint64_t> &values)
    {
        if(values.size()!=elements(shape))throw std::invalid_argument("canonical archive: integer dimensions");
        Bytes data;data.reserve(8*values.size());for(uint64_t x:values)appendU64(data,x);
        add(std::move(name),Type::U64,std::move(shape),std::move(data));
    }
    void reals(std::string name,std::vector<uint64_t> shape,const std::vector<double> &values)
    {
        if(values.size()!=elements(shape))throw std::invalid_argument("canonical archive: real dimensions");
        Bytes data;data.reserve(8*values.size());
        for(double x:values)
        {
            if(!std::isfinite(x))throw std::invalid_argument("canonical archive: nonfinite mechanical scalar");
            appendU64(data,doubleBits(x));
        }
        add(std::move(name),Type::F64,std::move(shape),std::move(data));
    }
    Bytes finish()const
    {
        Bytes out(magic.begin(),magic.end());appendU64(out,1);appendU64(out,generation_);appendU64(out,fields_.size());
        // Canonical lexicographic order, independent of registration order.
        for(const auto &[name,field]:fields_)
        {
            appendU64(out,name.size());out.insert(out.end(),name.begin(),name.end());
            out.push_back(uint8_t(field.type));out.push_back(uint8_t(field.dimensions.size()));
            for(uint64_t n:field.dimensions)appendU64(out,n);
            appendU64(out,field.data.size());out.insert(out.end(),field.data.begin(),field.data.end());
            if(out.size()>maximum_bytes-40)throw std::invalid_argument("canonical archive: payload byte bound");
        }
        const auto digest=sha256(out,out.size());out.insert(out.end(),digest.begin(),digest.end());
        out.insert(out.end(),footer.begin(),footer.end());return out;
    }
};

class Cursor
{
    const Bytes &bytes_;size_t limit_,offset_=0;
  public:
    Cursor(const Bytes &bytes,size_t limit):bytes_(bytes),limit_(limit)
    {if(limit>bytes.size())throw std::invalid_argument("canonical archive: cursor limit");}
    size_t offset()const{return offset_;}
    Bytes take(size_t size)
    {
        if(size>limit_-offset_)throw std::invalid_argument("canonical archive: truncated field");
        Bytes out(bytes_.begin()+offset_,bytes_.begin()+offset_+size);offset_+=size;return out;
    }
    uint8_t byte()
    {if(offset_==limit_)throw std::invalid_argument("canonical archive: truncated byte");return bytes_[offset_++];}
    uint64_t integer()
    {uint64_t out=0;for(size_t k=0;k<8;++k)out|=uint64_t(byte())<<(8*k);return out;}
};

/** A decoded candidate owns its bytes and never binds live simulation arrays.
 * Domain semantics/source/model/output checks belong to the later session
 * adapter. A valid byte envelope alone is deliberately NOT a restart pass.
 */
class Candidate
{
    Bytes canonical_;uint64_t generation_;
    std::map<std::string,Field> fields_;
    const Field &required(const std::string &name,Type type,const std::vector<uint64_t> &shape)const
    {
        const auto found=fields_.find(name);
        if(found==fields_.end() || found->second.type!=type || found->second.dimensions!=shape)
            throw std::invalid_argument("canonical archive: required field type/dimensions");
        return found->second;
    }
  public:
    const std::vector<uint64_t> &dimensions(const std::string &name)const
    {
        const auto found=fields_.find(name);
        if(found==fields_.end())throw std::invalid_argument("canonical archive: missing field dimensions");
        return found->second.dimensions;
    }
    Candidate(Bytes bytes,uint64_t explicit_generation):canonical_(std::move(bytes)),generation_(explicit_generation)
    {
        if(!generation_ || canonical_.size()<72 || canonical_.size()>maximum_bytes)
            throw std::invalid_argument("canonical archive: envelope size/generation");
        const size_t prefix=canonical_.size()-40;
        if(!std::equal(magic.begin(),magic.end(),canonical_.begin()) ||
           !std::equal(footer.begin(),footer.end(),canonical_.end()-8))
            throw std::invalid_argument("canonical archive: magic/footer");
        const auto digest=sha256(canonical_,prefix);
        if(!std::equal(digest.begin(),digest.end(),canonical_.begin()+prefix))
            throw std::invalid_argument("canonical archive: SHA256 mismatch");
        Cursor cursor(canonical_,prefix);cursor.take(8);
        if(cursor.integer()!=1 || cursor.integer()!=generation_)
            throw std::invalid_argument("canonical archive: schema or explicit generation mismatch");
        const uint64_t count=cursor.integer();
        if(count>maximum_fields)throw std::invalid_argument("canonical archive: field count");
        std::string previous;
        for(uint64_t i=0;i<count;++i)
        {
            const uint64_t length=cursor.integer();
            if(length==0 || length>maximum_name)throw std::invalid_argument("canonical archive: field name size");
            const auto encoded=cursor.take(length);std::string name(encoded.begin(),encoded.end());validName(name);
            if(!previous.empty() && name<=previous)throw std::invalid_argument("canonical archive: noncanonical field order");
            previous=name;
            const Type type=Type(cursor.byte());const auto rank=cursor.byte();width(type);
            if(rank>3)throw std::invalid_argument("canonical archive: rank");
            std::vector<uint64_t> shape;for(size_t k=0;k<rank;++k)shape.push_back(cursor.integer());
            const uint64_t size=cursor.integer();
            if(size!=elements(shape)*width(type))throw std::invalid_argument("canonical archive: field size");
            auto data=cursor.take(size);
            if(type==Type::F64)
            {
                Cursor values(data,data.size());
                for(uint64_t k=0;k<size/8;++k)if(!std::isfinite(fromBits(values.integer())))
                    throw std::invalid_argument("canonical archive: nonfinite mechanical scalar");
            }
            fields_.emplace(std::move(name),Field{type,std::move(shape),std::move(data)});
        }
        if(cursor.offset()!=prefix)throw std::invalid_argument("canonical archive: trailing payload");
    }
    void requireExactNames(const std::set<std::string> &expected)const
    {
        if(expected.size()!=fields_.size())throw std::invalid_argument("canonical archive: exact schema names");
        for(const auto &[name,field]:fields_)if(!expected.count(name))
            throw std::invalid_argument("canonical archive: unknown schema field");
    }
    const Bytes &canonicalBytes()const{return canonical_;}
    Bytes bytes(const std::string &name,const std::vector<uint64_t> &shape)const
    {return required(name,Type::Byte,shape).data;}
    std::vector<uint64_t> integers(const std::string &name,const std::vector<uint64_t> &shape)const
    {
        const auto &field=required(name,Type::U64,shape);Cursor cursor(field.data,field.data.size());
        std::vector<uint64_t> result;result.reserve(field.data.size()/8);
        while(cursor.offset()<field.data.size())result.push_back(cursor.integer());
        return result;
    }
    std::vector<double> reals(const std::string &name,const std::vector<uint64_t> &shape)const
    {
        const auto &field=required(name,Type::F64,shape);Cursor cursor(field.data,field.data.size());
        std::vector<double> result;result.reserve(field.data.size()/8);
        while(cursor.offset()<field.data.size())result.push_back(fromBits(cursor.integer()));
        return result;
    }
};
} // namespace SPH::peridynamics::canonical_archive
#endif
