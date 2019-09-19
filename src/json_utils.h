#ifndef json_utils_H_
#define json_utils_H_

#include <json_spirit.h>
#include <sstream>
#include <fstream>
#include <string>


#define VERBOSE 0

struct MembCarrier;
class Json_utils{
    public:
        const std::streampos getSize(const std::string&);
        const json_spirit::mValue& get_object_item(const json_spirit::mValue&, const std::string&);
        const json_spirit::mValue& get_array_item(const json_spirit::mValue&, size_t);
        MembCarrier JSONLoading(const std::string& filename);

    Json_utils(){};
    ~Json_utils(){};

};

struct MembCarrier{
    double* db_ptr, *db_ptr2; // db_ptr2 is meant to contain the k-space vector for 2D case!
    int* int_ptr;
    bool* boo_ptr;
    explicit MembCarrier(double*, double*, int*, bool*);
};

#endif /* json_utils_H_ */