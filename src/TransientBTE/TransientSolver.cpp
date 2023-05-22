//
// Created by yuehu on 2023/5/21.
//

#include "TransientBTE/transient.h"
#include <algorithm>
#include <iomanip>
#include <chrono>
using namespace std;

void Transient::solve(int Use_Backup, double error_temp_limit, double error_flux_limit,double deltaT, double totalT)
{
    double Num_Max_Iter=totalT/deltaT;
    _set_cell_matrix_larger();

    errorIncreaseTime=0;
    _set_initial(Use_Backup);
    
    
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++)
    {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local)
        {
            _get_bound_ee(iband_local,inf_local);
        }
    }
    _set_bound_ee_1();

    for (int ib = 0; ib < numBound; ++ib) {
        for (int inf = 0; inf < numDirectionLocal; ++inf) {
            for (int iband = 0; iband < numBandLocal; ++iband) {
                ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib*2]=temperature[0] * heatCapacity[matter[boundaryCell[ib][0]]][iband][inf];
                ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib*2+1]=temperature[0] * heatCapacity[matter[boundaryCell[ib][1]]][iband][inf];
            }
        }
    }
    for (int iband = 0; iband < numBand; ++iband) {
        for (int inf = 0; inf < numDirection ; ++inf) {
            for (int ib = 0; ib < numBound*2; ++ib) {
               // cout<< ebound[iband * numDirection * numBound * 2 + inf * numBound * 2 + ib]<<endl;
            }

        }

    }

    auto total_iter_time = chrono::microseconds(0);
    auto get_gradient_time = chrono::microseconds(0);
    auto get_Re_time=chrono::microseconds(0);
    auto solver1_time = chrono::microseconds(0);
    auto Boundary_time = chrono::microseconds(0);
    auto set_vertex_time=chrono::microseconds(0);
    auto face_time=chrono::microseconds(0);
    auto non_frourier_time = chrono::microseconds(0);
    auto set_bound_time = chrono::microseconds(0);
    auto macro_bound_time = chrono::microseconds(0);
    auto macro_iter_time = chrono::microseconds(0);
    auto trasfer1_time = chrono::microseconds(0);

    _get_CellMatrix_larger();
    ofstream outputT("TTG.dat");
    for (int nt = 0; nt < Num_Max_Iter; ++nt)
    {
        total_iter_time = chrono::microseconds(0);
        get_gradient_time = chrono::microseconds(0);
        get_Re_time=chrono::microseconds(0);
        solver1_time = chrono::microseconds(0);
        Boundary_time = chrono::microseconds(0);
        set_vertex_time=chrono::microseconds(0);
        face_time=chrono::microseconds(0);
        non_frourier_time = chrono::microseconds(0);
        set_bound_time = chrono::microseconds(0);
        macro_iter_time = chrono::microseconds(0);
        trasfer1_time = chrono::microseconds(0);

        auto total_iter_start = chrono::high_resolution_clock::now();
        copy();
        for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++)
        {
            for (int iband_local = 0; iband_local < numBandLocal; ++iband_local)
            {
                auto get_gradient_start =chrono::high_resolution_clock::now();
                _get_gradient_larger(0,iband_local,inf_local);
                auto get_gradient_end=chrono::high_resolution_clock::now();
                get_gradient_time+=chrono::duration_cast<chrono::microseconds>(get_gradient_end - get_gradient_start);

                auto get_Re_start =chrono::high_resolution_clock::now();
                _get_explicit_Re(nt,2,0,iband_local,inf_local,deltaT);
                auto get_Re_end =chrono::high_resolution_clock::now();
                get_Re_time+=chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                auto solve_start = chrono::high_resolution_clock::now();

                for (int icell = 0; icell < numCell; ++icell)
                {
                    //if(Re[icell]!=0)
                    //cout<<Re[icell]<<endl;
                    int inf = ((inf_local)*numProc + worldRank) % numDirection;
                    int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
                    energyDensity[iband_local][inf_local][icell] =energyDensity[iband_local][inf_local][icell] *(1 - deltaT / relaxationTime[matter[icell]][iband][inf]) - deltaT * Re[icell];
                }
                auto solve_end = chrono::high_resolution_clock::now();
                solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                auto boundary_start = chrono::high_resolution_clock::now();
                _get_bound_ee(iband_local,inf_local);
                auto boundary_end = chrono::high_resolution_clock::now();
                Boundary_time+=chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                 _recover_temperature(iband_local,inf_local);
                _get_total_energy(iband_local,inf_local);
                _get_heat_flux(iband_local,inf_local);


               MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            }
        }

        auto set_bound_start=chrono::high_resolution_clock::now();
        _set_bound_ee_1();
        auto set_bound_end=chrono::high_resolution_clock::now();
        set_bound_time+=chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);

        outputT<<(nt+1)*deltaT<<" "<<temperature[0]<<endl;
        auto total_iter_end = chrono::high_resolution_clock::now();
        total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
        MPI_Barrier(MPI_COMM_WORLD);

        if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit))
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
        }
        if (errorIncreaseTime >= 10000)
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
            if (worldRank == 0)
                cout << "error increases for 10 times, maybe the solution can not converge. Try \"-1\" in limiter" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (nt % 1 == 0)
        {
            if (worldRank == 0)
                _print_out();
        }

        if (worldRank == 0){

            cout << "  Time taken by inner loop: " << 1.0 * total_iter_time.count() / 1000 << " milliseconds" << endl;
            cout << "  Time taken by gradient 1: " << 1.0 * get_gradient_time.count() / 1000 << " milliseconds" << endl;
            cout << "  Time taken by BTE Re: " << 1.0 * get_Re_time.count() / 1000 << " milliseconds" << endl;
            cout << "  Time taken by BTE solver: " << 1.0 * solver1_time.count() /1000  << " milliseconds" << endl;
            cout << "  Time taken by Boundary: " << 1.0 * Boundary_time.count() / 1000 << " milliseconds"
            << " " << 1.0 * set_bound_time.count() / 1000 << " milliseconds" << endl;
            cout << "  Time taken by set_vertex: " << 1.0 * set_vertex_time.count() / 1000 << " milliseconds" << endl;
            cout << "  Time taken by face " << 1.0 * face_time.count() / 1000 << " milliseconds" << endl;
            cout << "  Time taken by non_fourier: " << 1.0 * non_frourier_time.count() / 1000 << " milliseconds" << endl;
            cout << "  Time taken by macro: " << 1.0 * macro_iter_time.count() / 1000 << " milliseconds" << endl;
            cout << "  Time taken by transfer: " << 1.0 * trasfer1_time.count() / 1000 << " milliseconds"<< endl;
            cout << "----------------------------------------------------------------------------------" << endl;
        }

        if(worldRank == 0)
        {
            for (int j = 1; j < numProc ; ++j) {
                MPI_Send(&nt,1,MPI_INT,j,10,MPI_COMM_WORLD);
            }

        }
        else
        {
            MPI_Status status;
            MPI_Recv(&nt,1,MPI_INT,0,10,MPI_COMM_WORLD,&status);
        }
        //cout<<nt<<endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
outputT.close();
    MPI_Barrier(MPI_COMM_WORLD);


    //_get_bound_temp();
    //_get_bound_flux();


#ifdef USE_TIME
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // cout << "end_iter1" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if (worldRank == 0)
        _print_out();
    if (worldRank == 0)

#ifdef USE_TIME
        cout << "Time taken by iteration: " << duration.count() * 0.001 << " milliseconds" << endl;
#endif

    _delete_cell_matrix();


}

void Transient::solve_first_order(int Use_Backup, int Use_Limiter, double error_temp_limit,
                                  double error_flux_limit, double deltaT, double totalT) {

}


