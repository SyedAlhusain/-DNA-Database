/********************************************
** File:    mytest.cpp
** Project: CMSC 341 Project 4, Spring 2022
** Author:  Syed Husain
** Date:    4/7/22
** E-mail:  ax18210@umbc.edu
** Desc:  This file test the dnadb object made in dnadb.cpp
** Course/Section : CMSC 341 
**/

#include "dnadb.cpp"
#include <random>
#include <vector>

enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL};
class Random {
public:
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor 
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else{ //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
    }
    void setSeed(int seedNum){
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum(){
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result*100.0)/100.0;
        return result;
    }
    
    private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};
class Tester{

    public:

    //test for testing all the test cases of dnadb object 
    bool DnadbAllTest();

    //test cases for insert()
    bool InsertNormalCase();
    bool InsertEdgeCase();
    bool InsertErrorCase();

    //test cases for getDNA()
    bool FindErrorCase();
    bool FindNormalCase();
    bool FindEdgeCaseNoRehash();

    //test cases for remove()
    bool RemoveNormalCase();
    bool RemoveEdgeCase();
    bool RemoveErrorCase();

    //test cases for rehash()
    bool RehashInsertNormal();
    bool RehashRemoveNormal();
    bool RehashInsertEdge();
    bool RehashRemoveEdge();

};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

int main(){
    // Preconditions: None
    // Postconditions: Call the Dnadb all test function test all the cases and then displays the result if 
    //the cases were passed or not

    Tester tester;

    //call the all test object pring program was successful else print it failed
    if(tester.DnadbAllTest() == true)
        cout << "\nThe Program worked successfully :) All the tests were passed!!" << endl;
    else
        cout << "\nThe Program failed :( All the tests did not passed succesfully !!" << endl;

}

unsigned int hashCode(const string str) {
   unsigned int val = 0 ;
   const unsigned int thirtyThree = 33 ;  // magic number from textbook
   for ( unsigned int i = 0 ; i < unsigned(str.length()); i++)
      val = val * thirtyThree + str[i] ;
   return val;
}

string sequencer(int size, int seedNum){
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0,3);
    rndObject.setSeed(seedNum);
    for (int i=0;i<size;i++){
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}


bool Tester::DnadbAllTest(){
    // Preconditions: None
    // Postconditions: Tests all the function to see if the program runs successfully

    Tester tester;

    //print the test and displays it result if the test fails return false and end the program

    cout << "\nTesting insert()(...), Normal case:\n\n";
    if (tester.InsertNormalCase() == true){
        cout << "\tNormal case of insert() passed!\n";
    }else{
        cout << "\tNormal case of insert() failed!\n"; 
        return false;
    }
    
    cout << "\nTesting insert()(...), Edge case:\n\n";
    if (tester.InsertEdgeCase() == true){
        cout << "\tEdge case of insert() passed!\n";
    }else{
        cout << "\tEdge case of insert() failed!\n"; 
        return false;
    }

     cout << "\nTesting insert()(...), Error case:\n\n";
    if (tester.InsertErrorCase() == true){
        cout << "\tError case of insert() passed!\n";
    }else{
        cout << "\tError case of insert() failed!\n"; 
        return false;
    }

    cout << "\nTesting find operation getDNA()(...) Error case:\n\n";
    if (tester.FindErrorCase() == true){
        cout << "\tError case case of getDNA() passed!\n";
    }else{
        cout << "\tError case case of getDNA() failed!\n"; 
        return false;
    }

    cout << "\nTesting find operation getDNA()(...) Normal case:\n\n";
    if (tester.FindNormalCase() == true){
        cout << "\tNormal case of getDNA() passed!\n";
    }else{
        cout << "\tNormal case of getDNA() failed!\n"; 
        return false;
    }

    cout << "\nTesting find operation getDNA()(...) Edge case No Rehash:\n\n";
    if (tester.FindEdgeCaseNoRehash() == true){
        cout << "\tEdge case case of getDNA() passed!\n";
    }else{
        cout << "\tEdge case case of getDNA() failed!\n"; 
        return false;
    }

    cout << "\nTesting remove()(...) Normal Case:\n\n";
    if (tester.RemoveNormalCase() == true){
        cout << "\tNormal Case of remove() passed!\n";
    }else{
        cout << "\tNormal Case of remove() failed!\n"; 
        return false;
    }
    
    cout << "\nTesting remove()(...) Edge Case:\n\n";
    if (tester.RemoveEdgeCase() == true){
        cout << "\tEdge Case of remove() passed!\n";
    }else{
        cout << "\tEdge Case of remove() failed!\n"; 
        return false;
    }

    cout << "\nTesting remove()(...) Error Case:\n\n";
    if (tester.RemoveErrorCase() == true){
        cout << "\tError Case of remove() passed!\n";
    }else{
        cout << "\tError Case of remove() failed!\n"; 
        return false;
    }

    cout << "\nTesting rehash of insert()(...) Normal Case:\n\n";
    if (tester.RehashInsertNormal() == true){
        cout << "\tNormal Case of rehashing insert() passed!\n";
    }else{
        cout << "\tNormal Case of rehashing insert()  failed!\n"; 
        return false;
    }

    cout << "\nTesting rehash of insert()(...) Edge Case:\n\n";
    if (tester.RehashInsertEdge() == true){
        cout << "\tEdge Case of rehashing insert() passed!\n";
    }else{
        cout << "\tEdge Case of rehashing insert()  failed!\n"; 
        return false;
    }

    cout << "\nTesting rehash of remove()(...) Normal Case:\n\n";
    if (tester.RehashRemoveNormal() == true){
        cout << "\tNormal Case of rehashing remove() passed!\n";
    }else{
        cout << "\tNormal Case of rehashing remove()  failed!\n"; 
        return false;
    }

    cout << "\nTesting rehash of remove()(...) Edge Case:\n\n";
    if (tester.RehashRemoveEdge() == true){
        cout << "\tEdge Case of rehashing remove() passed!\n";
    }else{
        cout << "\tEdge Case of rehashing remove()  failed!\n"; 
        return false;
    }
    
    return true; //return true if none of the test failed

}


bool Tester::InsertNormalCase(){
    // Preconditions: None
    // Postconditions: Return true if multiple non colliding key were inserted

    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    //iterates over and insert dna objects
    for (int i=0;i<49;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());

        //computes hash index to check if inserted correctly
        int hash_ind = (dnadb.m_hash)(dataObj.m_sequence);    
        int index = hash_ind % signed(dnadb.m_currentCap);

        //insert 
        if( i < 1){
            dnadb.insert(dataObj);
            dataList.push_back(dataObj);
        }else{
            //check for duplicate key and only insert non colliding key
            if(dnadb.checkKey(dataObj.m_sequence) != true){
                //insert object and check if return false
                dataList.push_back(dataObj);
                if(dnadb.insert(dataObj) == false){
                    return false;
                }
            }
        }

        //check if all node were inserted in correct index
        if(dnadb.m_currentTable[index].m_sequence == ""){
            return false;
        }
    }
    
    //check if all data was inserted
    bool result = true;
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++){
        result = result && (*it == dnadb.getDNA((*it).getSequence(), (*it).getLocId()));
    }

    //return false it some data are missing
    if (!result){
        return false;
    }

    //return false if colliding keys were inserted
    if(dnadb.m_currentSize == 49){
        return false;
    }

    return true; //return true once it passes all the test

}


bool Tester::InsertEdgeCase(){
    // Preconditions: None
    // Postconditions: Return true if colliding keys were inserted

    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    //iterates over and insert dna objects
    for (int i=0;i<20;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        //computes hash index to check if inserted correctly
        int hash_ind = (dnadb.m_hash)(dataObj.m_sequence);    
        int index = hash_ind % signed(dnadb.m_currentCap);

        //insert object and check if return false
        dataList.push_back(dataObj);
        
        if(dnadb.insert(dataObj) == false){
            return false;
        }

        //insert the colliding/duplicate key with different location
        dataObj.m_location =  RndLocation.getRandNum();
        //insert object and check if return false
        dataList.push_back(dataObj);
        
        if(dnadb.insert(dataObj) == false){
            return false;
        }

         //check if all node were inserted in correct index
        if(dnadb.m_currentTable[index].m_sequence == ""){
            return false;
        }

    }

    cout << "Line 370" << endl;
    dnadb.dump();

    //check if all data was inserted
    bool result = true;
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++){
        result = result && (*it == dnadb.getDNA((*it).getSequence(), (*it).getLocId()));
    }

    //return false it some data are missing
    if (!result){
        return false;
    }

    //check if any colliding key was not inserted
    if(dnadb.m_currentSize != 40){
        return false;
    }

    return true; //return true once it passes all the test

}

bool Tester::InsertErrorCase(){
    // Preconditions: None
    // Postconditions: Return true if duplicate dna is not inserted

    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    DNA duplicate_temp;

    //iterates over and insert dna objects
    for (int i=0;i<49;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        //computes hash index to check if inserted correctly
        int hash_ind = (dnadb.m_hash)(dataObj.m_sequence);    
        int index = hash_ind % signed(dnadb.m_currentCap);

        //insert object and check if return false
        if(dnadb.insert(dataObj) == false){
            return false;
        }

        //set the temp to an existing object
        if(i == 25){
            duplicate_temp = dataObj;
        }

        //check if all node were inserted in correct index
        if(dnadb.m_currentTable[index].m_sequence == ""){
            return false;
        }

    }

    //tries to insert duplicate object and return false if it inserted true
    if(dnadb.insert(duplicate_temp) == true){
        return false;
    }

    //check if the size wasn't increase by the insertion either
    if(dnadb.m_currentSize != 49){
        return false;
    }

    return true; //return true once it passes all the test

}

bool Tester::FindErrorCase(){
    // Preconditions: None
    // Postconditions: Return true, if DNA object does not exist in the database.

    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i=0;i<49;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        //insert object and check if return false        
        if(dnadb.insert(dataObj) == false){
            return false;
        }
    }

    //set a random temp
    DNA dataObj = DNA(sequencer(5, 50), RndLocation.getRandNum());

    //check if random temp that wasn't inserted exist in the table
    if(dnadb.getDNA(dataObj.m_sequence, dataObj.m_location).m_sequence != ""){
        //return false if it does, else return false
        return false;
    }

    return true; //return true once it passes all the test

}


bool Tester::FindNormalCase(){
    // Preconditions: None
    // Postconditions: Returns true if a DNA is found with a table consisting of non colliding keys
    
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    DNA temp;

    for (int i=0;i<20;i++){

        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());

        //insert the firt dna object
        if( i < 1){
            dnadb.insert(dataObj);
        }else{

            //check for duplicate key and only insert non colliding key
            if(dnadb.checkKey(dataObj.m_sequence) != true){
                temp = dataObj;
                //insert object and check if return false        
                if(dnadb.insert(dataObj) == false){
                    return false;
                }
            }
        }

    }

    //set the dna object to find the temp dna inserted 
    DNA dnaget = dnadb.getDNA(temp.getSequence(),temp.getLocId());

    //check if it return the correct temp object
    if(temp.getLocId() != dnaget.getLocId() && temp.getSequence() != dnaget.getSequence()){
        return false;
    }

    //check if the object is not empty
    if(dnadb.getDNA(dnaget.m_sequence, dnaget.m_location).m_sequence == ""){
        return false;
    }

    return true; //return true once it passes all the test

}
bool Tester::FindEdgeCaseNoRehash(){
    // Preconditions: None
    // Postconditions: Return true if DNA object if found with table consisting of multiple colliding keys

    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    DNA temp;

    for (int i=0;i<20;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());

        if(i == 7){
            temp = dataObj;

        }
        
        if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
            return false;
        }

        //insert the colliding/duplicate key with different location
        dataObj.m_location =  RndLocation.getRandNum();
        dnadb.insert(dataObj);

    }

    //set the dna object to find the temp dna inserted 
    DNA dnafind = dnadb.getDNA(temp.m_sequence, temp.m_location);

    //check if it return the correct temp object
    if(dnafind.m_sequence != temp.m_sequence && dnafind.m_location != temp.m_location){
        return false;
    }

    //check if the object is not empty
    if(dnadb.getDNA(dnafind.m_sequence, dnafind.m_location).m_sequence == ""){
        return false;
    }

    return true; //return true once it passes all the test

}

bool Tester::RemoveNormalCase(){
    // Preconditions: None
    // Postconditions: Returns true if a dna object is removed with a table consisting few non-colliding keys else return false

    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    DNA temp;

    for (int i=0;i<49;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());

        //insert the first object
        if( i < 1){   
            if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
                return false;}
            dataList.push_back(dataObj);
        }else{
             //check for duplicate key and only insert non colliding key
            if(dnadb.checkKey(dataObj.m_sequence) != true){
                temp = dataObj;
                if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
                    return false;}
                dataList.push_back(dataObj);
            }
        }
    }

    //check if temp object is exist in the table
    if(dnadb.getDNA(temp.m_sequence, temp.m_location).m_sequence == ""){
        return false;
    }

    //remove the temp
    if(dnadb.remove(temp) == false){
        return false;}

    //check if temp object is removed succesfully
    if(dnadb.getDNA(temp.m_sequence, temp.m_location).m_sequence != ""){
        return false;
    }

    bool result = true;
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++){
        result = result && (*it == dnadb.getDNA((*it).getSequence(), (*it).getLocId()));
    }

    //check if the temp object is removed and the data is missing
    if (result){
        return false;
    }

    return true; //return true once it passes all the test
}


bool Tester::RemoveEdgeCase(){
    // Preconditions: None
    // Postconditions: Returns true if a dna object is removed with a table consisting of colliding keys 
    // without being triggered else return false

    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    DnaDb dnadb_Duplicates(MINPRIME, hashCode);
    DNA temp;

    for (int i=0;i<8;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());

        dataList.push_back(dataObj);
        if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
            return false;}

        //insert the colliding/duplicate key with different location
        dataObj.m_location =  RndLocation.getRandNum();
        dataList.push_back(dataObj);
        if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
            return false;}
        
        //set the temp to be removed
        if(i == 3){
            temp = dataObj;
        }

    }

    //check if temp object is exist in the table
    if(dnadb.getDNA(temp.m_sequence, temp.m_location).m_sequence == ""){
        return false;
    }

    //remove the temp
    if(dnadb.remove(temp) == false){
        return false;}

        //check if temp object is removed succesfully
    if(dnadb.getDNA(temp.m_sequence, temp.m_location).m_sequence != ""){
        return false;
    }


    bool result = true;
    for (vector<DNA>::iterator it = dataList.begin(); it != dataList.end(); it++){
        result = result && (*it == dnadb.getDNA((*it).getSequence(), (*it).getLocId()));
    }

    //check if the temp object is removed and the data is missing
    if (result){
        return false;
    }

    return true; //return true once it passes all the test

}

bool Tester::RemoveErrorCase(){
    // Preconditions: None
    // Postconditions: Return true if non existing object is not removed(remove return false) or else return false

    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    DNA non_existing_temp;

    //iterates over and insert dna objects
    for (int i=0;i<49;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        //computes hash index to check if inserted correctly
        int hash_ind = (dnadb.m_hash)(dataObj.m_sequence);    
        int index = hash_ind % signed(dnadb.m_currentCap);

        //insert object and check if return false
        if(dnadb.insert(dataObj) == false){
            return false;
        }

        //check if all node were inserted in correct index
        if(dnadb.m_currentTable[index].m_sequence == ""){
            return false;
        }

    }

    non_existing_temp = DNA(sequencer(5, 50), RndLocation.getRandNum());

    //tries to remove a non existing object and return false if it removed true
    if(dnadb.remove(non_existing_temp) == true){
        return false;
    }

    //check if the size wasn't increase by the insertion either
    if(dnadb.m_currentSize != 49){
        return false;
    }

    return true; //return true once it passes all the test

}

bool Tester::RehashInsertNormal(){
    // Preconditions: None
    // Postconditions: Returns true if rehashing is triggered after a descent number of data insertion else returns false

    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    DNA temp;

    bool result = false;

    for (int i=0;i<52;i++){

        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
            return false;}

        //check if loading factor occurs closer to 0.5
        if(dnadb.lambda() > 0.49){
            result = true;
        }

        //check if anypoint new table is created reducing the size
        if(signed(dnadb.m_currentSize) != (i + 1)){
            result = true;
        }
    }

    //check if old table is set for rehashing tiggered
    if(dnadb.m_oldTable == nullptr){
        return false;
    }

    return result; //return true once it passes all the test
}


bool Tester::RehashInsertEdge(){
    // Preconditions: None
    // Postconditions: Returns true if after triggering rehash due to load factor, all live data 
    //is transferred to the new table and the old table is removed else return false


    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    DNA temp;

    bool result = false;

    for (int i=0;i<54;i++){

        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
            return false;}

        //check if loading factor occurs closer to 0.5
        if(dnadb.lambda() > 0.49){
            result = true;
        }

        //check if anypoint new table is created reducing the size
        if(signed(dnadb.m_currentSize) != (i + 1)){
            result = true;
        }

    }

    //check if old table is set for rehashing tiggered
    if(dnadb.m_oldTable == nullptr){
        return false;
    }

    DNA dataObj = DNA(sequencer(5, 55), RndLocation.getRandNum());
    if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
        return false;}

    //check if old table is deallocated after rehashing is tiggered
    if(dnadb.m_oldTable != nullptr){
        result = false;
    }

    return result; //return true once it passes all the test

}

bool Tester::RehashRemoveNormal() {
    // Preconditions: None
    // Postconditions: Return true if deleted ratio triggers a reshash after a descent number of data removal 
    //else returns false

    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i=0;i<50;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        dataList.push_back(dataObj);  // saving data for later use
        if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
            return false;}
    }

    //remove 80% of the node triggering a rehash
    for(int i= 0; i <  41; i++){
        if(dnadb.remove(dataList[i]) == false){ //return false if any removal return false
            return false;
        }
    }

    //check if a new table is allocated after triggering reshash
    if(dnadb.m_oldTable == nullptr){
        return false;
    }
    
    return true; //return true once it passes all the test
}

bool Tester::RehashRemoveEdge(){
    // Preconditions: None
    // Postconditions: Returns true if after triggering rehash due to delete ratio, i.e. 
    //all live data is transferred to the new table and the old table is removed else returns false

    vector<DNA> dataList;
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);

    for (int i=0;i<50;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        dataList.push_back(dataObj);  // saving data for later use
        if(dnadb.insert(dataObj) == false){ //insert object and check if return false        
            return false;}
    }

    //remove more than 80% of the node for triggering and completing the rehash
    for(int i= 0; i <  42; i++){
        if(dnadb.remove(dataList[i]) == false){ //return false if any removal return false
            return false;
        }
    }

    //check if the old table is deallocated after rehashing
    if(dnadb.m_oldTable != nullptr){
        return false;
    }

    return true; //return true once it passes all the test
    
}
