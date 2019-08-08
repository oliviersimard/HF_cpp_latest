#include "json_utils.h"


const std::streampos Json_utils::getSize(const std::string& filename){

    std::streampos begin, end;
    std::ifstream myfile(filename, std::ifstream::in);
    begin = myfile.tellg();
    myfile.seekg(0, std::ifstream::end);
    end = myfile.tellg();
    myfile.close();

    return end-begin;
}

const json_spirit::mValue& Json_utils::get_object_item(const json_spirit::mValue& element, const std::string& name){

    return element.get_obj().at(name);
}

const json_spirit::mValue& Json_utils::get_array_item(const json_spirit::mValue& element, size_t index){

    return element.get_array().at(index);
}

MembCarrier::MembCarrier(double* dp, double* dp2, int* ip) : db_ptr(dp), db_ptr2(dp2), int_ptr(ip){

}

MembCarrier Json_utils::JSONLoading(std::string& filename){

    std::ifstream JSONText(filename, std::ifstream::in);
    std::stringstream buffer;
    buffer << JSONText.rdbuf();
    json_spirit::mValue value;
    json_spirit::read(buffer,value); // Setting mValue object from buffer content.
    const std::string fileText(buffer.str()); // needs c_str() if printed using stdout.
    JSONText.close();
    if (VERBOSE > 0) std::cout << fileText.c_str() << std::endl;

    // Reading from json file
    //double precision
    const auto& Umax_val = get_object_item(value, "Umax");
    const auto& Ustep_val = get_object_item(value, "Ustep");
    const auto& Umin_val = get_object_item(value, "Umin");
    const auto& Vmax_val = get_object_item(value, "Vmax");
    const auto& Vstep_val = get_object_item(value, "Vstep");
    const auto& Vmin_val = get_object_item(value, "Vmin");
    const auto& betamax_val = get_object_item(value, "betamax");
    const auto& betastep_val = get_object_item(value, "betastep");
    const auto& betamin_val = get_object_item(value, "betamin");
    const auto& q_1D_val = get_object_item(value, "q_1D");
    //integers
    const auto& Nomega_val = get_object_item(value, "Nomega");
    const auto& N_it_val = get_object_item(value, "N_it");
    const auto& dims_val = get_object_item(value, "dims");
    const auto& gridK_val = get_object_item(value, "gridK");
    //array
    const auto& array_val = get_object_item(value, "q_2D");
    double container[2]; // Way to extract array inputs in json file!
    for (size_t i=0; i < array_val.get_array().size(); i++){
        container[i] = array_val.get_array().at(i).get_real();
    }

    if (VERBOSE > 0){
        std::cout << "Size of input file: " << getSize("params.json") << "\n";
    }
    
    //Collecting the json variables to instantiate param object
    double Umax = Umax_val.get_real(); double Vmax = Vmax_val.get_real();
    double Ustep = Ustep_val.get_real(); double Vstep = Vstep_val.get_real();
    double Umin = Umin_val.get_real(); double Vmin = Vmin_val.get_real();
    double betamax = betamax_val.get_real(); double betastep = betastep_val.get_real();
    double betamin = betamin_val.get_real();
    double q_1D = q_1D_val.get_real();
    int N_it = N_it_val.get_int(); int Nomega = Nomega_val.get_int();
    int dims = dims_val.get_int();
    int gridK = gridK_val.get_int();

    //Creating the arrays
    double dub[10] = { Umax, Ustep, Umin, Vmax, Vstep, Vmin, betamax, betastep, betamin, q_1D };
    // for (size_t i=0; i<=sizeof(dub)/sizeof(double); i++){
    //     std::cout << dub[i] << std::endl;
    // }
    int integ[4] = { Nomega, N_it, dims, gridK };
    // for (size_t i=0; i<=sizeof(integ)/sizeof(int); i++){
    //     std::cout << integ[i] << std::endl;
    // }

    return MembCarrier(dub,container,integ);
}