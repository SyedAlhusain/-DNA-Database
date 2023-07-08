/********************************************
** File:    dnadb.cpp
** Project: CMSC 341 Project 4, Spring 2022
** Author:  Syed Husain
** Date:    4/7/22
** E-mail:  ax18210@umbc.edu
** Desc:  This file makes the dnadb object to run by mytest.cpp
** Course/Section : CMSC 341 
**/

#include "dnadb.h"

DnaDb::DnaDb(int size, hash_fn hash){
    // Preconditions: None
    // Postconditions: Initializes variable and set the table

    //sets the hash
    m_hash = hash;   

    //check if size it ouside bound and set it to min or max depending on it value
    if(size > MAXPRIME)
        m_currentCap = MAXPRIME;
    else if(size < MINPRIME)
        m_currentCap = MINPRIME;
    else
    //if it doesn't go outside bound set it to current cap
        m_currentCap = size;

    //initalizes current table variables
    m_currentTable = new DNA [unsigned(m_currentCap)];
    m_currentSize = 0;              
    m_currNumDeleted = 0; 

    //initalizes old table variables
    m_oldTable = nullptr;   
    m_oldCap = 0;  
    m_oldSize = 0;                 
    m_oldNumDeleted = 0;  
    
}

DnaDb::~DnaDb(){
    // Preconditions: None
    // Postconditions: The old and current table are both deallocated by deleting all the nodes 

    //check if current table exist
    if(m_currentTable != nullptr){

        //deletes the table and reinitializes the variables
        delete[] m_currentTable;
        m_currentTable = nullptr;
        m_currentCap = 0;

    }

    //check if old table exist
    if(m_oldTable != nullptr){

        //deletes the table and reinitializes the variables
        delete[] m_oldTable ;
        m_oldTable = nullptr;
        m_oldCap = 0;     
    } 
    

}

bool DnaDb::checkDna(DNA dna){
    // Preconditions: None
    // Postconditions: Return true if a dna already exist in the table else return true

    //iterates over the table check for dna with exact same sequence and location
    if (m_currentTable != nullptr){
        for (int i = 0; i < signed(m_currentCap); i++) {
            //once found returns true else returns false
            if(m_currentTable[i].m_sequence == dna.m_sequence && m_currentTable[i].m_location == dna.m_location)
                return true;
        }
    }

    return false;

}


bool DnaDb::insert(DNA dna){
    // Preconditions: None
    // Postconditions: Return true if it inserts an object into the current hash table else returns false

    //check if m_location goes outside bounds and returns false
    if(dna.m_location > MAXLOCID || dna.m_location < MINLOCID)
        return false;

    //check if its a duplicate dna and return false
    if(checkDna(dna) == true){
        return false;
    }

    //bool to return at the end if an object is succesfully inserted or not
    bool found = false;

    //find the index to insert
    int hash_ind = (m_hash)(dna.m_sequence);    
    int index = hash_ind % signed(m_currentCap);

    //check if the node is empty or deleted
    if(m_currentTable[index].m_sequence == "" || m_currentTable[index].m_sequence == "DELETED"){

        //if its an deleted node decrement the counter since it already filled
        if(m_currentTable[index].m_sequence == "DELETED"){
            m_currentSize -= 1;
            m_currNumDeleted -= 1;
        }

        //insert the dna to it hashed index and increment the size
        m_currentTable[index] = dna;
        found = true;
        m_currentSize += 1;

    //if an hash index is already filled looks for the next hash index
    }else{
        bool flag = true;
        int temp_index;
        int i = 0;

        //iter over the table and find the next hashed index
        while(flag == true && i < signed(m_currentCap)){

            //use quadratic probing formula to find the next hashed index
            temp_index = (index + (i * i)) % signed(m_currentCap);

            //check if the node is empty or deleted
            if(m_currentTable[temp_index].m_sequence == "" || m_currentTable[temp_index].m_sequence == "DELETED"){

                //if its an deleted node decrement the curr since it already filled and curr deleted size since it active
                if(m_currentTable[temp_index].m_sequence == "DELETED"){
                    m_currentSize -= 1;
                    m_currNumDeleted -= 1;
                }

                //insert the dna to it hashed index and increment the size
                m_currentTable[temp_index] = dna;

                //set the found to true to return at the end and set flag to false to end the loop
                found = true;
                flag = false;

            }
            i++;
        }

        m_currentSize += 1;
    }

    //Reshashes if old table exist and keeps on reshaing until its deallocated
    if(m_oldTable != nullptr){
        RehashTable();
    }

    //Start the reshaing process once load or deleted factor exceed their trigger value
    if(lambda() > 0.5f){
        RehashTable();
    }

    //check if all node in the old table are deleted and added to current table 
    if(m_oldNumDeleted == m_oldSize){

        //once it does check if the table still exist and dellocates the old table
        if (m_oldTable != nullptr){
            delete[] m_oldTable ;

            //Reinitializes the variables
            m_oldTable = nullptr;
            m_oldCap = 0;
            m_oldSize = 0;
            m_oldNumDeleted = 0;
        }
    } 

    //return true or false depending on whether node was inserted or no
    return found;

}

void DnaDb::RehashTable(){
    // Preconditions: Load or deleted factor must exceed their trigger value
    // Postconditions: Transfer all the node 25% incrementally each time

    //if there no old table to transfer create a new table and swap it memories
    if(m_oldCap < 1 ){

        //equation to determine the capcity of new table
        int tempcap = (4 * (m_currentSize - m_currNumDeleted));
        //set the table capacity to next prime number by passing the equation
        m_oldCap = findNextPrime(tempcap);

        //allocates the new table and swaps it memories
        m_oldTable = new DNA [unsigned(m_oldCap)];
        swap(m_currentTable, m_oldTable);
        swap(m_currentCap, m_oldCap);
        swap(m_oldSize, m_currentSize);
        swap(m_currNumDeleted, m_oldNumDeleted);

    }

    int del_25 = 0; //variable to count 25% of the node
    int i = 0; //variable to count the whole table

    //checks if old table exist
    if (m_oldTable != nullptr){

        //loops until it transfer 25% of the node or reaches the end of current table
        while( del_25 < signed(m_oldSize / 4) && i < signed(m_oldCap)){
            
            //find the hash index to insert
            int hash_ind2 = (m_hash)(m_oldTable[i].m_sequence);    
            int index2 = hash_ind2 % signed(m_currentCap);

            //check if there an empty space to insert
            if(m_currentTable[index2].m_sequence == ""){

                //check if the node inserted it not deleted or empty
                if(m_oldTable[i].m_location != 0 && m_oldTable[i].m_sequence != ""){

                    //insert the node from the old table(after marking it deleted) into the current table 
                    m_currentTable[index2] = m_oldTable[i];
                    m_oldTable[i] = DELETED;

                    //increment the old deleted size and curr size along with del 25 to count 25% node
                    m_oldNumDeleted++;
                    m_currentSize++;
                    del_25++;

                }
            
            }else{

                //variables
                bool flag = true; //flag to end the loop
                int temp_index2; //hash index
                int i2 = 1; //count the table value

                //iterates over the loop and inserts the node or find the end of table
                while(flag == true && i2 < signed(m_currentCap)){

                    temp_index2 = (index2 + (i2 * i2)) % signed(m_currentCap); //hash index to insert

                    //check if there an empty space to insert
                    if(m_currentTable[temp_index2].m_sequence == ""){     

                        //check if the node inserted it not deleted or empty
                        if(m_oldTable[i].m_location != 0 && m_oldTable[i].m_sequence != ""){ 

                            //insert the node from the old table(after marking it deleted) into the current table                   
                            m_currentTable[temp_index2] = m_oldTable[i];
                            m_oldTable[i] = DELETED;

                            //increment the old deleted size and curr size along with del 25 to count 25% node
                            m_oldNumDeleted++;
                            del_25++;
                            m_currentSize++;
                            flag = false;
                            
                        }
                        
                    }
                    i2++;
                }
            }
            i++; 
        }   
    }

}


bool DnaDb::remove(DNA dna){

    // Preconditions: Current or old Table should exist
    // Postconditions: Return true if removes a data point from either the current hash table or 
    //the old hash table where the object is stored else returns false

    //bool to return at the end if an object is succesfully removed or not
    bool found = false;

    //check if current table exist
    if (m_currentTable != nullptr){

        //hash index to find the index to search on
        int hash_ind = (m_hash)(dna.m_sequence);    
        int index = hash_ind % signed(m_currentCap);

        bool flag = true; //bool to end the loop
        int temp_index; //index of the dna
        int i = 0; //to increment loopnig through the table

        //loops until it find the dna or reaches the end of table
        while(flag == true && i < signed(m_currentCap)){
            //index to search
            temp_index = (index + (i * i)) % signed(m_currentCap);

            //checks if the node have same locationa dn sequence with the dna to be removed
            if(m_currentTable[temp_index].m_location == dna.m_location && m_currentTable[temp_index].m_sequence == dna.m_sequence){

                //increment the curr delete and set it to deleted
                m_currNumDeleted += 1;
                m_currentTable[temp_index] = DELETED;

                //set found to true and flag to false in order to end the loop
                found = true;
                flag = false;
            }
            i++;
        }

    }

    //check if old table exist
    if (m_oldTable != nullptr){

        //hash index to find the index to search on
        int hash_ind = (m_hash)(dna.m_sequence);    
        int index = hash_ind % signed(m_oldCap);

        bool flag = true; //bool to end the loop
        int temp_index; //index of the dna
        int i = 0; //to increment loopnig through the table

        //loops until it find the dna or reaches the end of table
        while(flag == true && i < signed(m_oldCap)){
            //index to search
            temp_index = (index + (i * i)) % signed(m_oldCap);

            //checks if the node have same locationa dn sequence with the dna to be removed
            if(m_oldTable[temp_index].m_location == dna.m_location && m_oldTable[temp_index].m_sequence == dna.m_sequence){

                //increment the old delete and set it to deleted
                m_oldNumDeleted += 1;
                m_oldTable[temp_index] = DELETED;

                //set found to true and flag to false in order to end the loop
                found = true;
                flag = false;
            }
            i++;
        }
    }

    //check if all node in the old table are deleted and added to current table 
    if(m_oldNumDeleted == m_oldSize){

        //once it does check if the table still exist and dellocates the old table
        if (m_oldTable != nullptr){
            delete[] m_oldTable ;

            //Reinitializes the variables
            m_oldTable = nullptr;
            m_oldCap = 0;
            m_oldSize = 0;
            m_oldNumDeleted = 0;
        }
    } 

    //Reshashes if old table exist and keeps on reshaing until its deallocated
    if(m_oldTable != nullptr){
        RehashTable();
    }

    //Start the reshaing process once load or deleted factor exceed their trigger value
    if(deletedRatio() > 0.8f){
        RehashTable();
    }

    return found; //return true or false depending on whether node was removed or not

}


DNA DnaDb::getDNA(string sequence, int location){   
    // Preconditions: Current or old Table should exist
    // Postconditions: Returns dna object found else returns an empty object

    //create an empty object
    DNA temp("", 0);

    //check if current table exist
    if (m_currentTable != nullptr){

        //hash index to find the index to search on
        int hash_ind = (m_hash)(sequence);    
        int index = hash_ind % signed(m_currentCap);

        bool flag = true; //bool to end the loop
        int temp_index; //index of the dna
        int i = 0; //to increment loopnig through the table

        //loops until it find the dna or reaches the end of table
        while(flag == true && i < signed(m_currentCap)){
            //index to search
            temp_index = (index + (i * i)) % signed(m_currentCap);

             //if both dna matches set the empty object to found node
            if(m_currentTable[temp_index].m_sequence == sequence && m_currentTable[temp_index].m_location == location){
                temp = m_currentTable[temp_index];
                //set found to true and flag to false in order to end the loop
                flag = false;
            }
            i++;
        }

    }

    //check if old table exist
    if (m_oldTable != nullptr){

        //hash index to find the index to search on
        int hash_ind = (m_hash)(sequence);    
        int index = hash_ind % signed(m_oldCap);

        bool flag = true; //bool to end the loop
        int temp_index; //index of the dna
        int i = 0; //to increment loopnig through the table

        //loops until it find the dna or reaches the end of table
        while(flag == true && i < signed(m_oldCap)){
            //index to search
            temp_index = (index + (i * i)) % signed(m_oldCap);

            //if both dna matches set the empty object to found node
            if(m_oldTable[temp_index].m_location == location && m_oldTable[temp_index].m_sequence == sequence){
                temp = m_oldTable[temp_index];
                //set found to true and flag to false in order to end the loop
                flag = false;
            }
            i++;
        }
    }

    //return the dna object if found else returns an empty object
    return temp;
    
}

float DnaDb::lambda() const {  
    // Preconditions: None
    // Postconditions: Returns the load factor of the current hash table. 

    float ratio; //float variable to return

    //calculates the ratio of occupied buckets to the table capacity. 
    ratio = float(m_currentSize) / float(m_currentCap);

    return ratio; //return the ratio
  
}

float DnaDb::deletedRatio() const {
    // Preconditions: None
    // Postconditions: Returns the deleted ratio factor of the current hash table. 

    float ratio; //float variable to return

    //calculates the deleted buckets to the total number of occupied buckets .
    ratio = float(m_currNumDeleted) / float(m_currentSize);

    return ratio; //return the ratio
    
}

void DnaDb::dump() const {
    // Preconditions: None
    // Postconditions: dumps the contents of the current hash table and the old has table if it exists. 

    cout << "Dump for current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < signed(m_currentCap); i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }

    cout << "Dump for old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < signed(m_oldCap); i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}


bool DnaDb::isPrime(int number){
    // Preconditions: None
    // Postconditions: Check if the number passed is prime
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current){
    // Preconditions: Current table must exist
    // Postconditions: Find the next prime

    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) { 
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0) 
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

bool DnaDb::checkKey(string key){
    // Preconditions: None
    // Postconditions: Return true if it a duplicate key else return false

    //iterates over the table and check if the key matches a dna sequence 
    if (m_currentTable != nullptr){
        for (int i = 0; i < signed(m_currentCap); i++) {
            //return ture if it matches else return false
            if(m_currentTable[i].m_sequence == key)
                return true;
                
        }
    }
    return false;

}

DNA::DNA(string sequence, int location) {
    if ((location >= MINLOCID && location <= MAXLOCID) ||
        (location == 0 && sequence == "DELETED")){
        // this is a normal or a DELETED object
        m_sequence = sequence;
        m_location = location;
    }
    else{
        // this is the empty object
        m_sequence = "";
        m_location = 0;
    }
}

string DNA::getSequence() const {
    return m_sequence;
}

int DNA::getLocId() const {
    return m_location;
}

// Overloaded assignment operator
const DNA& DNA::operator=(const DNA& rhs){
    if (this != &rhs){
        m_sequence = rhs.m_sequence;
        m_location = rhs.m_location;
    }
    return *this;
}

// Overloaded insertion operator.  Prints DNA's sequence (key),
// and the location ID. This is a friend function in DNA class.
ostream& operator<<(ostream& sout, const DNA &dna ) {
    if (!dna.m_sequence.empty())
        sout << dna.m_sequence << " (Location ID " << dna.m_location << ")";
    else
        sout << "";
    return sout;
}

// Overloaded equality operator. This is a friend function in DNA class.
// To test inequality we may negate the results of this operator.
bool operator==(const DNA& lhs, const DNA& rhs){
    return ((lhs.m_sequence == rhs.m_sequence) && (lhs.m_location == rhs.m_location));

}
