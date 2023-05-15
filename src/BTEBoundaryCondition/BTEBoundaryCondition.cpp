//
// Created by yuehu on 2021/9/8.
//

#include "BTEBoundaryCondition/BTEBoundaryCondition.h"
#include "BTEBand/BTEBand.h"
#include "sstream"
using namespace std;
BTEBoundaryCondition::BTEBoundaryCondition(ifstream& inFile,BTEBand *bands) {
    if (!inFile.is_open()) {
        cout << "DEfBUG: BC file not open" << endl;
        exit(1);
    }
    int num;
    string line;
    inFile >> num;
    for (int i = 0; i < 3; i++)
        getline(inFile, line);
    for (int i = 0; i < num; i++) {
        BoundaryCondition bc;
        int temp;
        inFile >> bc.index  >> bc.type;
        if (bc.type==4||bc.type==5)
        {
            inFile>>bc.temperature>>bc.Name_Interface_file;
        } else
        {
            inFile>> bc.temperature;
        }

        getline(inFile, line);
        boundaryConditions.push_back(bc);
    }
    for (int i = 0; i < num ; ++i) {
        if (boundaryConditions[i].type==4||boundaryConditions[i].type==5)
        {
            ifstream Ininterface(boundaryConditions[i].Name_Interface_file);
            if (!Ininterface.is_open()) {
                cout << "DEBUG: interface file not open" << endl;
                exit(1);
            }
            while(getline(Ininterface, line)){
                //cout<<line.find("matter")<<" "<<line.length()<<endl;
                if(line.find("matter")>0 &&line.find("matter")< line.length()){
                    //cout<<line<<endl;
                    transcoeff trans;
                    trans.matter=stoi(line.substr(0, line.find(" ")));

                    while (getline(Ininterface, line)) {
                        //cout<<line<<endl;
                        if (line.length() <= 1) {
                            break;
                        }
                        stringstream ss(line);
                        double temp;
                        ss>>temp;
                        trans.trans.push_back(temp);
                    }
                    if (trans.trans.size()!=bands->bands.size())
                    {
                        cout<<"mismatch band vs interface"<<endl;
                    }
                    boundaryConditions[i].matter_trans.push_back(trans);
                }
            }
            Ininterface.close();
        }

    }
}
BTEBoundaryCondition::~BTEBoundaryCondition(){
    cout<<"~BTEBoundaryCondition is activated!!"<<endl;
}