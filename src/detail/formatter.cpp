// String formatter to make runtime_errors look good
// taken from: http://stackoverflow.com/questions/12261915/howto-throw-stdexceptions-with-variable-messages

#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      class formatter
      {
      public:
          formatter() {}
          ~formatter() {}
      
          template <typename Type>
          formatter & operator << (const Type & value)
          {
              stream_ << value;
              return *this;
          }
      
          std::string str() const         { return stream_.str(); }
          operator std::string () const   { return stream_.str(); }
      
          enum ConvertToString 
          {
              to_str
          };
          std::string operator >> (ConvertToString) { return stream_.str(); }
      
      private:
          std::stringstream stream_;
      
          formatter(const formatter &);
          formatter & operator = (formatter &);
      };
    };
  };
};
