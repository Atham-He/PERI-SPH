#ifndef PD_CASE_CONFIG_H
#define PD_CASE_CONFIG_H
#include "base_data_type.h"
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace SPH::peridynamics
{
/** Strict key=value case inputs. Framework command-line options are preserved. */
class CaseParameters
{
  public:
    void bind(const std::string &key, Real &value)
    {
        setters_[key] = [&value](const std::string &s)
        {
            size_t end; Real parsed = std::stod(s, &end);
            if (end != s.size() || !std::isfinite(parsed)) throw std::invalid_argument("expected finite number");
            value = parsed;
        };
        getters_[key] = [&value]() { return value; };
    }
    void bind(const std::string &key, int &value)
    {
        setters_[key] = [&value](const std::string &s)
        {
            size_t end; int parsed = std::stoi(s, &end);
            if (end != s.size()) throw std::invalid_argument("expected integer");
            value = parsed;
        };
        getters_[key] = [&value]() { return Real(value); };
    }
    void bind(const std::string &key, bool &value)
    {
        setters_[key] = [&value](const std::string &s)
        {
            if (s == "true" || s == "1") value = true;
            else if (s == "false" || s == "0") value = false;
            else throw std::invalid_argument("expected true or false");
        };
        getters_[key] = [&value]() { return Real(value); };
    }
    void parse(int &argc, char **argv)
    {
        std::string file; int output = 1;
        for (int i = 1; i < argc; ++i)
        {
            std::string option(argv[i]);
            if (option == "--case-config" || option.rfind("--case-config=", 0) == 0)
            {
                if (!file.empty()) throw std::invalid_argument("only one --case-config is accepted");
                if (option == "--case-config")
                {
                    if (++i == argc) throw std::invalid_argument("--case-config requires a path");
                    file = argv[i];
                }
                else file = option.substr(14);
            }
            else argv[output++] = argv[i];
        }
        argc = output; argv[output] = nullptr;
        if (file.empty()) return;
        std::ifstream input(file);
        if (!input) throw std::invalid_argument("cannot open case configuration: " + file);
        std::map<std::string, bool> seen;
        std::string line; size_t number = 0;
        while (std::getline(input, line))
        {
            ++number; line = trim(line.substr(0, line.find('#')));
            if (line.empty()) continue;
            size_t equal = line.find('=');
            if (equal == std::string::npos) throw std::invalid_argument("expected key=value at line " + std::to_string(number));
            std::string key = trim(line.substr(0, equal)), value = trim(line.substr(equal + 1));
            if (!setters_.count(key)) throw std::invalid_argument("unknown case parameter: " + key);
            if (seen[key]) throw std::invalid_argument("duplicate case parameter: " + key);
            seen[key] = true;
            try { setters_.at(key)(value); }
            catch (const std::exception &e) { throw std::invalid_argument(key + ": " + e.what()); }
        }
    }
    void save(const std::string &path) const
    {
        std::ofstream output(path);
        if (!output) throw std::runtime_error("cannot write effective parameters: " + path);
        output << std::setprecision(17);
        for (const auto &entry : getters_) output << entry.first << "=" << entry.second() << "\n";
    }

  private:
    static std::string trim(const std::string &s)
    {
        size_t first = s.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        return s.substr(first, s.find_last_not_of(" \t\r\n") - first + 1);
    }
    std::map<std::string, std::function<void(const std::string &)>> setters_;
    std::map<std::string, std::function<Real()>> getters_;
};
}
#endif
