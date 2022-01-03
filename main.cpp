#include<iostream>
#include<chrono>
#include<ctime>
#include<thread>
#include"include/ocm.h"
using namespace std;
void my_binary(int64_t n);
uint64_t cal(string str_k_mer);
uint64_t reverse_compliment(uint64_t cal_kmer, int kmer_length);
template <typename ccmbase_obj> void update_count_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name);
template <typename ccmbase_obj> void update_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round);
template <typename ccmbase_obj> void update_count_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round, int total_round);
//std::chrono::time_point<std::chrono::steady_clock> start_time,end_time;
float total_time = 0;

string sample_kmer="AAACCGGTCTGTTGGGACCACT";  // len = 22
string INPUT_FASTA_FILE = "",OUTPUT_SKETCH_FILE="", query_file_name="", query_result_file_name="output/query_output.csv";   //  "rymv.sim.fa"; //
unsigned int kmer_len = 0;
int counter=0;
unsigned int NP = 20, NH = 7, TOTAL_ROUND = 4,SEED=137, counter_w,NUM_THREAD;
bool CANONICALIZE = true, CONSERVATIVE = false;

int main(int argc, char *argv[]){
    string mode(argv[1]);
    if( mode == "test"){
        for(int i= 2; i<argc-1; i++){
            string arg(argv[i]);
            //string param(argv[++i]);
            if(arg=="-k") kmer_len=stoi(argv[++i]);
            else if(arg=="-h") NH =stoi(argv[++i]);
            else if(arg=="-w") NP = log2(stoi(argv[++i]));
            else if(arg=="-n") TOTAL_ROUND=stoi(argv[++i]);
            else if(arg=="-t") NUM_THREAD=stoi(argv[++i]);
            else if(arg=="-o") OUTPUT_SKETCH_FILE = argv[++i];
            else if(arg=="-q") query_file_name = argv[++i];
            else if(arg=="-fa") INPUT_FASTA_FILE=argv[++i];
            else if(arg=="-r") CANONICALIZE = false;
            else if(arg=="-c") CONSERVATIVE = true;
        }
        if(kmer_len <= 0 || INPUT_FASTA_FILE=="" || OUTPUT_SKETCH_FILE==""){
            cout<<"INVALID INPUT"; return 0;
        }
        //end inouting
        //cout<<NP<<" "<<NH<<" "<<TOTAL_ROUND<<" "<<OUTPUT_SKETCH_FILE<<" "<<INPUT_FASTA_FILE<<" "<<CANONICALIZE<<" ";
        if(!CONSERVATIVE){
            //construct OCM
            total_time = 0;
            sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash,2> sketch1(NP,NH,137);
            for(int r = 0; r< TOTAL_ROUND; r++){
                if (r > 0){
                    // for all kmers update collision
                    auto start_time = std::chrono::high_resolution_clock::now();
                    update_collision_from_file(INPUT_FASTA_FILE,kmer_len,sketch1,r);
                    auto end_time = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<float> duration = end_time - start_time;
                    total_time+=duration.count();
                    cout<<"Updating collision for ocms round "<<r<<" is completed in "<<duration.count()<<" seconds\n";
                }

                sketch1.clear_core();
                // for all kmer update count.
                auto start_time = std::chrono::high_resolution_clock::now();
                //sketch1.update_count_from_file(INPUT_FASTA_FILE,kmer_len,CANONICALIZE);
                /// updating count through thread

                std::vector<std::thread> threads;
                for (std::size_t i = 0; i < NUM_THREAD; i++) {
                    string ith_filename = INPUT_FASTA_FILE+"_"+to_string(i)+".fa";
                    cout<<"Trying to open "<<ith_filename<<endl;
                    threads.emplace_back(&sketch::ocm::ocmbase<uint64_t, sketch::hash::WangHash,2>::update_count_from_file, &sketch1, ith_filename,kmer_len,CANONICALIZE);
                }

                for (auto& thread : threads) {
                    thread.join();
                }
                /// thread ends

                auto end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<float> duration = end_time - start_time;
                total_time+=duration.count();
                cout<<"Updating count for ocms round "<<r<<" is completed in "<<duration.count()<<" seconds\n";
                //sketch1.showCounters(cal(sample_kmer));
            }
            cout<<"Constructing ocms is completed in :"<<total_time<<" seconds\n";
            //end construct OCM
            //sketch1.save_sketch(OUTPUT_SKETCH_FILE);


            //query starts
            ifstream infile(query_file_name);
            if(!infile.good()){cout<<"Couldn't open input query file\n";}
            string kmer;
            int true_count;
            ofstream query_result_file;
            query_result_file.open(query_result_file_name, ios::out);
            if(!query_result_file.good()){cout<<"Couldn't open output query file file\n";}
            query_result_file<<"kmer,true_count,estimated_count\n";
            while (infile >> kmer >> true_count){
                //cout<<kmer<<","<<true_count<<","<<query_sketch.est_count(cal(kmer))<<endl;
                query_result_file <<kmer<<","<<true_count<<","<<sketch1.est_count(cal(kmer))<<endl;
            }
            //query ends

        }
        else{
        //construct occm
        total_time = 0;
        sketch::ocm::ocmbase <uint64_t, sketch::hash::WangHash,2> sketch2(NP,NH,137);
        for(int r = 0; r< TOTAL_ROUND; r++){
            if (r > 0){
                // for all kmers update collision
                auto start_time = std::chrono::high_resolution_clock::now();
                update_collision_from_file(INPUT_FASTA_FILE,kmer_len,sketch2,r);
                auto end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<float> duration = end_time - start_time;
                total_time+=duration.count();
                cout<<"Updating collision for occms round "<<r<<" is completed in "<<duration.count()<<" seconds\n";
            }
            sketch2.clear_core();
            // for all kmer update count collision
            auto start_time = std::chrono::high_resolution_clock::now();
            //for(auto kmer: kmers_vec) sketch2.update_count_collision(kmer,r,TOTAL_ROUND);
            update_count_collision_from_file(INPUT_FASTA_FILE,kmer_len,sketch2,r,TOTAL_ROUND);
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> duration = end_time - start_time;
            total_time+=duration.count();
            cout<<"Updating count-collision for occms round "<<r<<" is completed in "<<duration.count()<<" seconds\n";
            //sketch2.showCounters(cal(sample_kmer));

        }
        //end occm
        cout<<"Constructing occm is completed in :"<<total_time<<" seconds\n";
        //sketch2.save_sketch(OUTPUT_SKETCH_FILE);

        //query starts
            ifstream infile(query_file_name);
            if(!infile.good()){cout<<"Couldn't open input query file\n";}
            string kmer;
            int true_count;
            ofstream query_result_file;
            query_result_file.open(query_result_file_name, ios::out);
            if(!query_result_file.good()){cout<<"Couldn't open output query file file\n";}
            query_result_file<<"kmer,true_count,estimated_count\n";
            while (infile >> kmer >> true_count){
                //cout<<kmer<<","<<true_count<<","<<query_sketch.est_count(cal(kmer))<<endl;
                query_result_file <<kmer<<","<<true_count<<","<<sketch2.est_count(cal(kmer))<<endl;
            }
            //query ends
        }
    }

    return 0;
}

void my_binary(int64_t n) //converting k-mer to binary representation
{
    int a[100];
    string str="";
    int i;
    for(i=0; i<64; i++)
    {
        if(n>0)
        {
            a[i]=n%2;
            n= n/2;
        }
        else a[i]=0;

    }
    cout<<endl<<"Binary: ";
    for(int j=i-1 ;j>=0 ;j--)
    {
        cout<<a[j];
        //str+=to_string(a[i]);
    }
    cout<<endl;
}




uint64_t cal(string str_k_mer)
{
    const char* char_array = str_k_mer.c_str();
    uint64_t k_mer = 0;
    for(int j=0; char_array[j]!='\0'; j++)
    {
        switch(char_array[j])
        {
            case 'A':
                    k_mer = k_mer<<2; // A=00
                    break;
            case 'T':
                    k_mer = k_mer<<2;
                    k_mer = k_mer | 1;   //T=01
                    break;
            case 'G':
                    k_mer = k_mer<<2;  //G=10
                    k_mer = k_mer | 2;
                    break;
            case 'C':
                    k_mer = k_mer<<2; //C=11
                    k_mer = k_mer | 3;
                    break;
        }
    }
    //cout<<"cal value: "<<k_mer;
    //my_binary(k_mer);
    return k_mer;

}

uint64_t reverse_compliment(uint64_t cal_kmer, int kmer_length)
{
    uint64_t k_mer = 0;
    uint64_t mask = 3;

    for(int i=0; i<kmer_length; i++)
    {
        switch(cal_kmer & mask)
        {

        case 0:
            k_mer = k_mer<<2;
            k_mer = k_mer | 1;   //A=00->T=01
            break;
        case 1:
            k_mer = k_mer<<2; //T=01->A=00
            break;
        case 2:
            k_mer = k_mer<<2; //G=10->C=11
            k_mer = k_mer | 3;
            break;
        case 3:
            k_mer = k_mer<<2;  //C=11->G=10
            k_mer = k_mer | 2;
            break;

        }
        cal_kmer=cal_kmer>>2;
    }
    return k_mer;
}




template <typename ccmbase_obj> void update_count_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name)
{
    std::ifstream input(file);
    if(!input.good()){cout<<"Can't open"<<endl;}

    std::string line, name, content;
    while( std::getline( input, line ).good()){ //reading 1 line from input

        if( line.empty() || line[0] == '>' ){
            if( !name.empty() ){
                int l=content.length();
                for(int i=0; i<=l-len_k_mer; i++) //extracting k-mers from the line
                {
                    string part = content.substr(i, len_k_mer);
                    int64_t k_mer = cal(part);
                    /////////Updating part.
                    ccmbase_obj_name.update_count(k_mer);
                    if(CANONICALIZE)ccmbase_obj_name.update_count(reverse_compliment(k_mer,kmer_len));
                    if(k_mer == cal(sample_kmer))counter++;
                }
                name.clear();
            }
            if( !line.empty() ){name = line.substr(1);}
            content.clear();
        }
        else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){
                name.clear();
                content.clear();
            }
            else {
                content += line;
            }
        }
    }
    if( !name.empty() ){
        //std::cout << name << " : " << content << std::endl;
        int l=content.length();
        for(int i=0; i<=l-len_k_mer; i++)
        {
            string part = content.substr(i, len_k_mer);
            int64_t k_mer = cal(part);

            //// updating part.
            ccmbase_obj_name.update_count(k_mer);
            if(CANONICALIZE)ccmbase_obj_name.update_count(reverse_compliment(k_mer,kmer_len));
            if(k_mer == cal(sample_kmer))counter++; //manual count of the sample k-mer


        }
    }
}


template <typename ccmbase_obj> void update_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round)
{
    std::ifstream input(file);
    if(!input.good()){cout<<"Can't open"<<endl;}

    std::string line, name, content;
    while( std::getline( input, line ).good()){ //reading 1 line from input

        if( line.empty() || line[0] == '>' ){
            if( !name.empty() ){
                int l=content.length();
                for(int i=0; i<=l-len_k_mer; i++) //extracting k-mers from the line
                {
                    string part = content.substr(i, len_k_mer);
                    int64_t k_mer = cal(part);
                    /////////Updating part.
                    ccmbase_obj_name.update_collision(k_mer, round);
                    if(CANONICALIZE)ccmbase_obj_name.update_collision(reverse_compliment(k_mer,kmer_len), round);
                }
                name.clear();
            }
            if( !line.empty() ){name = line.substr(1);}
            content.clear();
        }
        else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){
                name.clear();
                content.clear();
            }
            else {
                content += line;
            }
        }
    }
    if( !name.empty() ){
        //std::cout << name << " : " << content << std::endl;
        int l=content.length();
        for(int i=0; i<=l-len_k_mer; i++)
        {
            string part = content.substr(i, len_k_mer);
            int64_t k_mer = cal(part);

            //// updating part.
            ccmbase_obj_name.update_collision(k_mer, round);
            if(CANONICALIZE)ccmbase_obj_name.update_collision(reverse_compliment(k_mer,kmer_len), round);

        }
    }
}

template <typename ccmbase_obj> void update_count_collision_from_file(string file, int len_k_mer, ccmbase_obj &ccmbase_obj_name, int round, int total_round)
{
    std::ifstream input(file);
    if(!input.good()){cout<<"Can't open"<<endl;}

    std::string line, name, content;
    while( std::getline( input, line ).good()){ //reading 1 line from input

        if( line.empty() || line[0] == '>' ){
            if( !name.empty() ){
                int l=content.length();
                for(int i=0; i<=l-len_k_mer; i++) //extracting k-mers from the line
                {
                    string part = content.substr(i, len_k_mer);
                    int64_t k_mer = cal(part);
                    /////////Updating part.
                    ccmbase_obj_name.update_count_collision(k_mer, round, total_round);
                    if(CANONICALIZE)ccmbase_obj_name.update_count_collision(reverse_compliment(k_mer,kmer_len), round, total_round);
                }
                name.clear();
            }
            if( !line.empty() ){name = line.substr(1);}
            content.clear();
        }
        else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){
                name.clear();
                content.clear();
            }
            else {
                content += line;
            }
        }
    }
    if( !name.empty() ){
        //std::cout << name << " : " << content << std::endl;
        int l=content.length();
        for(int i=0; i<=l-len_k_mer; i++)
        {
            string part = content.substr(i, len_k_mer);
            int64_t k_mer = cal(part);
            //// updating part.
            ccmbase_obj_name.update_count_collision(k_mer, round, total_round);
            if(CANONICALIZE)ccmbase_obj_name.update_count_collision(reverse_compliment(k_mer,kmer_len), round, total_round);
        }
    }
}
