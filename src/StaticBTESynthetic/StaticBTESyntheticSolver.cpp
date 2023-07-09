//
// Created by yuehu on 2023/4/7.
//
#include "StaticBTESynthetic/StaticBTESynthetic.h"
#include <algorithm>
#include <iomanip>
#include <chrono>
using namespace std;



void StaticBTESynthetic::solve(int Use_Backup, int Num_Max_Iter, int Use_Limiter, double error_temp_limit, double error_flux_limit)
{

    _set_face_matrix();
    _set_cell_matrix_larger();
    _set_matrix("LU");
    errorIncreaseTime=0;
    _set_initial(Use_Backup);
    if (worldRank==0)
    {cout<<"Begin to calculate coefficient ..."<<endl;}
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            _get_coefficient(iband_local, inf_local);
        }
    }
    if (worldRank==0)
    {cout<<"Finish calculating coefficient"<<endl;}
    _get_coefficient_macro();


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

    _get_face_matrix();
    _get_CellMatrix_larger();
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            _get_gradient_larger(Use_Limiter,iband_local,inf_local);
            _get_bound_ee(iband_local,inf_local);
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                  numBound * 2,
                  MPI_DOUBLE,
                  (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                  numBound * 2,
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);
        }
    }
    _set_bound_ee_1();
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
                _get_gradient_larger(Use_Limiter,iband_local,inf_local);
                auto get_gradient_end=chrono::high_resolution_clock::now();
                get_gradient_time+=chrono::duration_cast<chrono::microseconds>(get_gradient_end - get_gradient_start);

                auto get_Re_start =chrono::high_resolution_clock::now();
                _get_Re(iband_local,inf_local);

                auto get_Re_end =chrono::high_resolution_clock::now();
                get_Re_time+=chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                auto solve_start = chrono::high_resolution_clock::now();

                sol1=cg[iband_local][inf_local].solve(Re1);
                for (int i = 0; i < numCell; ++i)
                {
                    energyDensity[iband_local][inf_local][i]= sol1(i) + energyDensity[iband_local][inf_local][i];
                }
                auto solve_end = chrono::high_resolution_clock::now();
                solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                auto boundary_start = chrono::high_resolution_clock::now();
                _get_bound_ee(iband_local,inf_local);
                auto boundary_end = chrono::high_resolution_clock::now();
                Boundary_time+=chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                auto vertex_start = chrono::high_resolution_clock::now();
                _set_vertex_energydensity(iband_local,inf_local);

                auto vertex_end = chrono::high_resolution_clock::now();
                set_vertex_time+=chrono::duration_cast<chrono::microseconds>(vertex_end - vertex_start);

                auto face_start = chrono::high_resolution_clock::now();
                _get_face_gradient(iband_local,inf_local);
                auto face_end = chrono::high_resolution_clock::now();
                face_time+=chrono::duration_cast<chrono::microseconds>(face_end - face_start);

                auto fourier_start = chrono::high_resolution_clock::now();
                _get_non_fourier(iband_local,inf_local);
                auto fourier_end = chrono::high_resolution_clock::now();
                non_frourier_time+=chrono::duration_cast<chrono::microseconds>(fourier_end - fourier_start);
                auto transfer_start = chrono::high_resolution_clock::now();
                _recover_temperature(iband_local,inf_local);
                _get_total_energy(iband_local,inf_local);
                _get_heat_flux(iband_local,inf_local);


                MPI_Allreduce(ReMacroLocal, ReMacro, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                auto transfer_end = chrono::high_resolution_clock::now();
                trasfer1_time+=chrono::duration_cast<chrono::microseconds>(transfer_end - transfer_start);
            }
        }

        auto set_bound_start=chrono::high_resolution_clock::now();
        _set_bound_ee_1();
        auto set_bound_end=chrono::high_resolution_clock::now();
        set_bound_time+=chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);

        auto set_macro_start=chrono::high_resolution_clock::now();
        _set_macro_boundary_condition();
        auto set_macro_end=chrono::high_resolution_clock::now();
        macro_bound_time+=chrono::duration_cast<chrono::microseconds>(set_macro_end-set_macro_start);

        auto macro_iter_start=chrono::high_resolution_clock::now();
        _marco_solution();
        auto macro_iter_end=chrono::high_resolution_clock::now();
        macro_iter_time+=chrono::duration_cast<chrono::microseconds>(macro_iter_end-macro_iter_start);

        auto total_iter_end = chrono::high_resolution_clock::now();
        total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
        MPI_Barrier(MPI_COMM_WORLD);

        if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit))
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
        }
        if (nt >= 5000 && nt != Num_Max_Iter && nt % 100 == 0)
        {
            if (worldRank == 0)
                cout << "More than 5,000 iterations have been performed. Please try \"Synthetic\" in IterativeScheme" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (errorIncreaseTime >= 10)
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
            if (worldRank == 0)
                cout << "Error increases for 10 times, maybe the solution can not converge. Try to set Order=1" << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);


    _get_bound_temp();
    _get_bound_flux();


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
    _delete_matrix("LU");
    _delete_face_matrix();
    _delete_cell_matrix();


}

void StaticBTESynthetic::solve_Iterative(int Use_Backup, int Num_Max_Iter, int Use_Limiter, double error_temp_limit, double error_flux_limit)
{

    _set_face_matrix();
    _set_cell_matrix_larger();
    _set_matrix("BICGSTAB");
    errorIncreaseTime=0;
    _set_initial(Use_Backup);

    _get_coefficient_macro();


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

    _get_face_matrix();
    _get_CellMatrix_larger();
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            _get_gradient_larger(Use_Limiter,iband_local,inf_local);
            _get_bound_ee(iband_local,inf_local);
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                          numBound * 2,
                          MPI_DOUBLE,
                          (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                          numBound * 2,
                          MPI_DOUBLE,
                          MPI_COMM_WORLD);
        }
    }
    _set_bound_ee_1();
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
                _get_gradient_larger(Use_Limiter,iband_local,inf_local);
                auto get_gradient_end=chrono::high_resolution_clock::now();
                get_gradient_time+=chrono::duration_cast<chrono::microseconds>(get_gradient_end - get_gradient_start);

                auto get_Re_start =chrono::high_resolution_clock::now();
                _get_Re(iband_local,inf_local);

                auto get_Re_end =chrono::high_resolution_clock::now();
                get_Re_time+=chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                auto solve_start = chrono::high_resolution_clock::now();

                if (worldRank==0 && inf_local==0 && iband_local==0) {
                    cout<<"Begin to calculate coefficient ..."<<endl;
                }
                _get_coefficient_Iterative(iband_local,inf_local);
                if (worldRank==0 && inf_local==numDirectionLocal-1 && iband_local==numBandLocal-1) {
                    cout<<"Finish calculating coefficient"<<endl;
                }

                Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double>> cg1;
                sol1=cg1.compute(stiffMatrix).solve(Re1);
                for (int i = 0; i < numCell; ++i)
                {
                    energyDensity[iband_local][inf_local][i]= sol1(i) + energyDensity[iband_local][inf_local][i];
                }
                auto solve_end = chrono::high_resolution_clock::now();
                solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                auto boundary_start = chrono::high_resolution_clock::now();
                _get_bound_ee(iband_local,inf_local);
                auto boundary_end = chrono::high_resolution_clock::now();
                Boundary_time+=chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                auto vertex_start = chrono::high_resolution_clock::now();
                _set_vertex_energydensity(iband_local,inf_local);

                auto vertex_end = chrono::high_resolution_clock::now();
                set_vertex_time+=chrono::duration_cast<chrono::microseconds>(vertex_end - vertex_start);

                auto face_start = chrono::high_resolution_clock::now();
                _get_face_gradient(iband_local,inf_local);
                auto face_end = chrono::high_resolution_clock::now();
                face_time+=chrono::duration_cast<chrono::microseconds>(face_end - face_start);

                auto fourier_start = chrono::high_resolution_clock::now();
                _get_non_fourier(iband_local,inf_local);
                auto fourier_end = chrono::high_resolution_clock::now();
                non_frourier_time+=chrono::duration_cast<chrono::microseconds>(fourier_end - fourier_start);
                auto transfer_start = chrono::high_resolution_clock::now();
                _recover_temperature(iband_local,inf_local);
                _get_total_energy(iband_local,inf_local);
                _get_heat_flux(iband_local,inf_local);


                MPI_Allreduce(ReMacroLocal, ReMacro, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                auto transfer_end = chrono::high_resolution_clock::now();
                trasfer1_time+=chrono::duration_cast<chrono::microseconds>(transfer_end - transfer_start);
            }
        }

        auto set_bound_start=chrono::high_resolution_clock::now();
        _set_bound_ee_1();
        auto set_bound_end=chrono::high_resolution_clock::now();
        set_bound_time+=chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);

        auto set_macro_start=chrono::high_resolution_clock::now();
        _set_macro_boundary_condition();
        auto set_macro_end=chrono::high_resolution_clock::now();
        macro_bound_time+=chrono::duration_cast<chrono::microseconds>(set_macro_end-set_macro_start);

        auto macro_iter_start=chrono::high_resolution_clock::now();
        _marco_solution();
        auto macro_iter_end=chrono::high_resolution_clock::now();
        macro_iter_time+=chrono::duration_cast<chrono::microseconds>(macro_iter_end-macro_iter_start);

        auto total_iter_end = chrono::high_resolution_clock::now();
        total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
        MPI_Barrier(MPI_COMM_WORLD);

        if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit))
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
        }
        if (nt >= 5000 && nt != Num_Max_Iter && nt % 100 == 0)
        {
            if (worldRank == 0)
                cout << "More than 5,000 iterations have been performed. Please try \"Synthetic\" in IterativeScheme" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (errorIncreaseTime >= 10)
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
            if (worldRank == 0)
                cout << "Error increases for 10 times, maybe the solution can not converge. Try to set Order=1" << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);


    _get_bound_temp();
    _get_bound_flux();


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
    _delete_matrix("BICGSTAB");
    _delete_face_matrix();
    _delete_cell_matrix();


}

void StaticBTESynthetic::solve_firstorder
(int Use_Backup, int Num_Max_Iter, int Use_Limiter,
 double error_temp_limit, double error_flux_limit)
 {


    _set_face_matrix();
    _set_matrix("LU");
    errorIncreaseTime=0;
    _set_initial(Use_Backup);
    if (worldRank==0)
    {cout<<"Begin to calculate coefficient ..."<<endl;}
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            _get_coefficient(iband_local, inf_local);
        }
    }
    if (worldRank==0)
    {cout<<"Finish calculating coefficient"<<endl;}
    _get_coefficient_macro();
    //_set_bound_ee_1();

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

    _get_face_matrix();
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            //_get_gradient_larger(Use_Limiter,iband_local,inf_local);
            _get_bound_ee_firstorder(iband_local,inf_local);
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                  numBound * 2,
                  MPI_DOUBLE,
                  (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                  numBound * 2,
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);
        }

    }
    _set_bound_ee_1();

     for (int i = 0; i < numCell; ++i) {
         gradientX[i]=0;
         gradientY[i]=0;
         gradientZ[i]=0;
         limit[i]=1;
     }

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
                auto get_Re_start =chrono::high_resolution_clock::now();
                _get_Re_firstorder(iband_local,inf_local);
                auto get_Re_end =chrono::high_resolution_clock::now();
                get_Re_time+=chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                auto solve_start = chrono::high_resolution_clock::now();

                sol1=cg[iband_local][inf_local].solve(Re1);
                for (int i = 0; i < numCell; ++i)
                {
                    energyDensity[iband_local][inf_local][i]= sol1(i);
                }
                auto solve_end = chrono::high_resolution_clock::now();
                solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                auto boundary_start = chrono::high_resolution_clock::now();
                _get_bound_ee_firstorder(iband_local,inf_local);
                auto boundary_end = chrono::high_resolution_clock::now();
                Boundary_time+=chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                auto vertex_start = chrono::high_resolution_clock::now();
                _set_vertex_energydensity(iband_local,inf_local);

                auto vertex_end = chrono::high_resolution_clock::now();
                set_vertex_time+=chrono::duration_cast<chrono::microseconds>(vertex_end - vertex_start);

                auto face_start = chrono::high_resolution_clock::now();
                _get_face_gradient(iband_local,inf_local);
                auto face_end = chrono::high_resolution_clock::now();
                face_time+=chrono::duration_cast<chrono::microseconds>(face_end - face_start);

                auto fourier_start = chrono::high_resolution_clock::now();
                _get_non_fourier(iband_local,inf_local);
                auto fourier_end = chrono::high_resolution_clock::now();
                non_frourier_time+=chrono::duration_cast<chrono::microseconds>(fourier_end - fourier_start);
                auto transfer_start = chrono::high_resolution_clock::now();
                _recover_temperature(iband_local,inf_local);
                _get_total_energy(iband_local,inf_local);
                _get_heat_flux(iband_local,inf_local);


                MPI_Allreduce(ReMacroLocal, ReMacro, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                auto transfer_end = chrono::high_resolution_clock::now();
                trasfer1_time+=chrono::duration_cast<chrono::microseconds>(transfer_end - transfer_start);
            }
        }

        auto set_bound_start=chrono::high_resolution_clock::now();
        _set_bound_ee_1();
        auto set_bound_end=chrono::high_resolution_clock::now();
        set_bound_time+=chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);

        auto set_macro_start=chrono::high_resolution_clock::now();
        _set_macro_boundary_condition();
        auto set_macro_end=chrono::high_resolution_clock::now();
        macro_bound_time+=chrono::duration_cast<chrono::microseconds>(set_macro_end-set_macro_start);

        auto macro_iter_start=chrono::high_resolution_clock::now();
        _marco_solution();
        auto macro_iter_end=chrono::high_resolution_clock::now();
        macro_iter_time+=chrono::duration_cast<chrono::microseconds>(macro_iter_end-macro_iter_start);

        auto total_iter_end = chrono::high_resolution_clock::now();
        total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
        MPI_Barrier(MPI_COMM_WORLD);

        if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit))
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
        }
        if (nt >= 5000 && nt != Num_Max_Iter && nt % 100 == 0)
        {
            if (worldRank == 0)
                cout << "More than 5,000 iterations have been performed. Please try \"Synthetic\" in IterativeScheme" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (errorIncreaseTime >= 10)
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
            if (worldRank == 0)
                cout << "Error increases for 10 times, maybe the solution can not converge. Try to set Order=1" << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);


    _get_bound_temp();
    _get_bound_flux();


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
    _delete_matrix("LU");
    _delete_face_matrix();


 }

 void StaticBTESynthetic::solve_firstorder_Iterative
 (int Use_Backup, int Num_Max_Iter, int Use_Limiter,
  double error_temp_limit, double error_flux_limit)
  {


    _set_face_matrix();
    _set_matrix("BICGSTAB");
    errorIncreaseTime=0;
    _set_initial(Use_Backup);

    _get_coefficient_macro();
    //_set_bound_ee_1();

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

    _get_face_matrix();
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            //_get_gradient_larger(Use_Limiter,iband_local,inf_local);
            _get_bound_ee_firstorder(iband_local,inf_local);
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                          numBound * 2,
                          MPI_DOUBLE,
                          (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                          numBound * 2,
                          MPI_DOUBLE,
                          MPI_COMM_WORLD);
        }

    }
    _set_bound_ee_1();

    for (int i = 0; i < numCell; ++i) {
        gradientX[i]=0;
        gradientY[i]=0;
        gradientZ[i]=0;
        limit[i]=1;
    }

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
                auto get_Re_start =chrono::high_resolution_clock::now();
                _get_Re_firstorder(iband_local,inf_local);
                auto get_Re_end =chrono::high_resolution_clock::now();
                get_Re_time+=chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                auto solve_start = chrono::high_resolution_clock::now();

                if (worldRank==0 && inf_local==0 && iband_local==0) {
                    cout<<"Begin to calculate coefficient ..."<<endl;
                }
                _get_coefficient_Iterative(iband_local,inf_local);
                if (worldRank==0 && inf_local==numDirectionLocal-1 && iband_local==numBandLocal-1) {
                    cout<<"Finish calculating coefficient"<<endl;
                }

                Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double>> cg1;
                sol1=cg1.compute(stiffMatrix).solve(Re1);
                for (int i = 0; i < numCell; ++i)
                {
                    energyDensity[iband_local][inf_local][i]= sol1(i) ;
                }
                auto solve_end = chrono::high_resolution_clock::now();
                solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                auto boundary_start = chrono::high_resolution_clock::now();
                _get_bound_ee_firstorder(iband_local,inf_local);
                auto boundary_end = chrono::high_resolution_clock::now();
                Boundary_time+=chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                auto vertex_start = chrono::high_resolution_clock::now();
                _set_vertex_energydensity(iband_local,inf_local);

                auto vertex_end = chrono::high_resolution_clock::now();
                set_vertex_time+=chrono::duration_cast<chrono::microseconds>(vertex_end - vertex_start);

                auto face_start = chrono::high_resolution_clock::now();
                _get_face_gradient(iband_local,inf_local);
                auto face_end = chrono::high_resolution_clock::now();
                face_time+=chrono::duration_cast<chrono::microseconds>(face_end - face_start);

                auto fourier_start = chrono::high_resolution_clock::now();
                _get_non_fourier(iband_local,inf_local);
                auto fourier_end = chrono::high_resolution_clock::now();
                non_frourier_time+=chrono::duration_cast<chrono::microseconds>(fourier_end - fourier_start);
                auto transfer_start = chrono::high_resolution_clock::now();
                _recover_temperature(iband_local,inf_local);
                _get_total_energy(iband_local,inf_local);
                _get_heat_flux(iband_local,inf_local);


                MPI_Allreduce(ReMacroLocal, ReMacro, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                auto transfer_end = chrono::high_resolution_clock::now();
                trasfer1_time+=chrono::duration_cast<chrono::microseconds>(transfer_end - transfer_start);
            }
        }

        auto set_bound_start=chrono::high_resolution_clock::now();
        _set_bound_ee_1();
        auto set_bound_end=chrono::high_resolution_clock::now();
        set_bound_time+=chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);

        auto set_macro_start=chrono::high_resolution_clock::now();
        _set_macro_boundary_condition();
        auto set_macro_end=chrono::high_resolution_clock::now();
        macro_bound_time+=chrono::duration_cast<chrono::microseconds>(set_macro_end-set_macro_start);

        auto macro_iter_start=chrono::high_resolution_clock::now();
        _marco_solution();
        auto macro_iter_end=chrono::high_resolution_clock::now();
        macro_iter_time+=chrono::duration_cast<chrono::microseconds>(macro_iter_end-macro_iter_start);

        auto total_iter_end = chrono::high_resolution_clock::now();
        total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
        MPI_Barrier(MPI_COMM_WORLD);

        if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit))
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
        }
        if (nt >= 5000 && nt != Num_Max_Iter && nt % 100 == 0)
        {
            if (worldRank == 0)
                cout << "More than 5,000 iterations have been performed. Please try \"Synthetic\" in IterativeScheme" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (errorIncreaseTime >= 10)
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
            if (worldRank == 0)
                cout << "Error increases for 10 times, maybe the solution can not converge. Try to set Order=1" << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);


    _get_bound_temp();
    _get_bound_flux();


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
    _delete_matrix("BICGSTAB");
    _delete_face_matrix();


  }

void StaticBTESynthetic::solve_DOM(int Use_Backup, int Num_Max_Iter, int Use_Limiter, double error_temp_limit, double error_flux_limit)
{



    _set_cell_matrix_larger();
    _set_matrix("LU");
    errorIncreaseTime=0;
    _set_initial(Use_Backup);
    if (worldRank==0)
    {cout<<"Begin to calculate coefficient ..."<<endl;}
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            _get_coefficient(iband_local, inf_local);
        }
    }
    if (worldRank==0)
    {cout<<"Finish calculating coefficient"<<endl;}
    //_set_bound_ee_1();

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
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            _get_gradient_larger(Use_Limiter,iband_local,inf_local);
            _get_bound_ee(iband_local,inf_local);
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                  numBound * 2,
                  MPI_DOUBLE,
                  (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                  numBound * 2,
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);
        }
    }
    _set_bound_ee_1();

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
                _get_gradient_larger(Use_Limiter,iband_local,inf_local);
                auto get_gradient_end=chrono::high_resolution_clock::now();
                get_gradient_time+=chrono::duration_cast<chrono::microseconds>(get_gradient_end - get_gradient_start);

                auto get_Re_start =chrono::high_resolution_clock::now();
                _get_Re(iband_local,inf_local);
                auto get_Re_end =chrono::high_resolution_clock::now();
                get_Re_time+=chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                auto solve_start = chrono::high_resolution_clock::now();

                sol1=cg[iband_local][inf_local].solve(Re1);
                for (int i = 0; i < numCell; ++i)
                {
                    energyDensity[iband_local][inf_local][i]= sol1(i) + energyDensity[iband_local][inf_local][i];
                }
                auto solve_end = chrono::high_resolution_clock::now();
                solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                auto boundary_start = chrono::high_resolution_clock::now();
                _get_bound_ee(iband_local,inf_local);
                auto boundary_end = chrono::high_resolution_clock::now();
                Boundary_time+=chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                 auto transfer_start = chrono::high_resolution_clock::now();
                _recover_temperature(iband_local,inf_local);
                _get_total_energy(iband_local,inf_local);
                _get_heat_flux(iband_local,inf_local);


                MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                auto transfer_end = chrono::high_resolution_clock::now();
                trasfer1_time+=chrono::duration_cast<chrono::microseconds>(transfer_end - transfer_start);
            }
        }

        auto set_bound_start=chrono::high_resolution_clock::now();
        _set_bound_ee_1();
        auto set_bound_end=chrono::high_resolution_clock::now();
        set_bound_time+=chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);


        auto total_iter_end = chrono::high_resolution_clock::now();
        total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
        MPI_Barrier(MPI_COMM_WORLD);

        if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit))
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
        }
        if (nt >= 5000 && nt != Num_Max_Iter && nt % 100 == 0)
        {
            if (worldRank == 0)
                cout << "More than 5,000 iterations have been performed. Please try \"Synthetic\" in IterativeScheme" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (errorIncreaseTime >= 10)
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
            if (worldRank == 0)
                cout << "Error increases for 10 times, maybe the solution can not converge. Try to set Order=1" << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);


    _get_bound_temp();
    _get_bound_flux();


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
    _delete_matrix("LU");

    _delete_cell_matrix();


}


void StaticBTESynthetic::solve_DOM_Iterative
(int Use_Backup, int Num_Max_Iter, int Use_Limiter,
 double error_temp_limit, double error_flux_limit)
 {



    _set_cell_matrix_larger();
    _set_matrix("BICGSTAB");
    errorIncreaseTime=0;
    _set_initial(Use_Backup);

    //_set_bound_ee_1();

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
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            _get_gradient_larger(Use_Limiter,iband_local,inf_local);
            _get_bound_ee(iband_local,inf_local);
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                          numBound * 2,
                          MPI_DOUBLE,
                          (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                          numBound * 2,
                          MPI_DOUBLE,
                          MPI_COMM_WORLD);
        }
    }
    _set_bound_ee_1();

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
                _get_gradient_larger(Use_Limiter,iband_local,inf_local);
                auto get_gradient_end=chrono::high_resolution_clock::now();
                get_gradient_time+=chrono::duration_cast<chrono::microseconds>(get_gradient_end - get_gradient_start);

                auto get_Re_start =chrono::high_resolution_clock::now();
                _get_Re(iband_local,inf_local);
                auto get_Re_end =chrono::high_resolution_clock::now();
                get_Re_time+=chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                auto solve_start = chrono::high_resolution_clock::now();

                if (worldRank==0 && inf_local==0 && iband_local==0) {
                    cout<<"Begin to calculate coefficient ..."<<endl;
                }
                _get_coefficient_Iterative(iband_local,inf_local);
                if (worldRank==0 && inf_local==numDirectionLocal-1 && iband_local==numBandLocal-1) {
                    cout<<"Finish calculating coefficient"<<endl;
                }

                Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double>> cg1;
                sol1=cg1.compute(stiffMatrix).solve(Re1);
                for (int i = 0; i < numCell; ++i)
                {
                    energyDensity[iband_local][inf_local][i]= sol1(i) + energyDensity[iband_local][inf_local][i];
                }
                auto solve_end = chrono::high_resolution_clock::now();
                solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                auto boundary_start = chrono::high_resolution_clock::now();
                _get_bound_ee(iband_local,inf_local);
                auto boundary_end = chrono::high_resolution_clock::now();
                Boundary_time+=chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                auto transfer_start = chrono::high_resolution_clock::now();
                _recover_temperature(iband_local,inf_local);
                _get_total_energy(iband_local,inf_local);
                _get_heat_flux(iband_local,inf_local);


                MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                auto transfer_end = chrono::high_resolution_clock::now();
                trasfer1_time+=chrono::duration_cast<chrono::microseconds>(transfer_end - transfer_start);
            }
        }

        auto set_bound_start=chrono::high_resolution_clock::now();
        _set_bound_ee_1();
        auto set_bound_end=chrono::high_resolution_clock::now();
        set_bound_time+=chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);


        auto total_iter_end = chrono::high_resolution_clock::now();
        total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
        MPI_Barrier(MPI_COMM_WORLD);

        if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit))
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
        }
        if (nt >= 5000 && nt != Num_Max_Iter && nt % 100 == 0)
        {
            if (worldRank == 0)
                cout << "More than 5,000 iterations have been performed. Please try \"Synthetic\" in IterativeScheme" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (errorIncreaseTime >= 10)
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
            if (worldRank == 0)
                cout << "Error increases for 10 times, maybe the solution can not converge. Try to set Order=1" << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);


    _get_bound_temp();
    _get_bound_flux();


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
    _delete_matrix("BICGSTAB");

    _delete_cell_matrix();


 }

void StaticBTESynthetic::solve_DOM_firstorder
(int Use_Backup, int Num_Max_Iter, int Use_Limiter,
double error_temp_limit, double error_flux_limit)
{


    _set_matrix("LU");
    errorIncreaseTime=0;
    _set_initial(Use_Backup);
    if (worldRank==0)
    {cout<<"Begin to calculate coefficient ..."<<endl;}
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            _get_coefficient(iband_local, inf_local);
        }
    }
    if (worldRank==0)
    {cout<<"Finish calculating coefficient"<<endl;}
    //_set_bound_ee_1();

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
    
    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            //_get_gradient_larger(Use_Limiter,iband_local,inf_local);
            _get_bound_ee_firstorder(iband_local,inf_local);
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                  numBound * 2,
                  MPI_DOUBLE,
                  (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                  numBound * 2,
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);
        }
    }
    _set_bound_ee_1();

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
                auto get_Re_start =chrono::high_resolution_clock::now();
                _get_Re_firstorder(iband_local,inf_local);
                auto get_Re_end =chrono::high_resolution_clock::now();
                get_Re_time+=chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                auto solve_start = chrono::high_resolution_clock::now();

                sol1=cg[iband_local][inf_local].solve(Re1);
                for (int i = 0; i < numCell; ++i)
                {
                    energyDensity[iband_local][inf_local][i]= sol1(i);
                }
                auto solve_end = chrono::high_resolution_clock::now();
                solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                auto boundary_start = chrono::high_resolution_clock::now();
                _get_bound_ee_firstorder(iband_local,inf_local);
                auto boundary_end = chrono::high_resolution_clock::now();
                Boundary_time+=chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                auto transfer_start = chrono::high_resolution_clock::now();
                _recover_temperature(iband_local,inf_local);
                _get_total_energy(iband_local,inf_local);
                _get_heat_flux(iband_local,inf_local);


                MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                auto transfer_end = chrono::high_resolution_clock::now();
                trasfer1_time+=chrono::duration_cast<chrono::microseconds>(transfer_end - transfer_start);
            }
        }

        auto set_bound_start=chrono::high_resolution_clock::now();
        _set_bound_ee_1();
        auto set_bound_end=chrono::high_resolution_clock::now();
        set_bound_time+=chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);


        auto total_iter_end = chrono::high_resolution_clock::now();
        total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
        MPI_Barrier(MPI_COMM_WORLD);

        if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit))
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
        }
        if (nt >= 5000 && nt != Num_Max_Iter && nt % 100 == 0)
        {
            if (worldRank == 0)
                cout << "More than 5,000 iterations have been performed. Please try \"Synthetic\" in IterativeScheme" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (errorIncreaseTime >= 10)
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
            if (worldRank == 0)
                cout << "Error increases for 10 times, maybe the solution can not converge. Try to set Order=1" << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);


    _get_bound_temp();
    _get_bound_flux();


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
    _delete_matrix("LU");



}

void StaticBTESynthetic::solve_DOM_firstorder_Iterative
(int Use_Backup, int Num_Max_Iter, int Use_Limiter,
 double error_temp_limit, double error_flux_limit)

 {


    _set_matrix("BICGSTAB");
    errorIncreaseTime=0;
    _set_initial(Use_Backup);

    //_set_bound_ee_1();

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

    for (int inf_local = 0; inf_local < numDirectionLocal; inf_local++) {
        for (int iband_local = 0; iband_local < numBandLocal; ++iband_local) {
            //_get_gradient_larger(Use_Limiter,iband_local,inf_local);
            _get_bound_ee_firstorder(iband_local,inf_local);
            int inf = ((inf_local) * numProc + worldRank) % numDirection;
            int iband = iband_local * (ceil(double(numProc) / double(numDirection))) + worldRank / numDirection;
            MPI_Allgather(eboundLocal + iband_local * numBound * 2,
                          numBound * 2,
                          MPI_DOUBLE,
                          (ebound + numBound * 2 * (inf - worldRank % numDirection)) + numDirection * numBound * 2 * (iband - worldRank / numDirection),
                          numBound * 2,
                          MPI_DOUBLE,
                          MPI_COMM_WORLD);
        }
    }
    _set_bound_ee_1();

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
                auto get_Re_start =chrono::high_resolution_clock::now();
                _get_Re_firstorder(iband_local,inf_local);
                auto get_Re_end =chrono::high_resolution_clock::now();
                get_Re_time+=chrono::duration_cast<chrono::microseconds>(get_Re_end - get_Re_start);

                auto solve_start = chrono::high_resolution_clock::now();

                if (worldRank==0 && inf_local==0 && iband_local==0) {
                    cout<<"Begin to calculate coefficient ..."<<endl;
                }
                _get_coefficient_Iterative(iband_local,inf_local);
                if (worldRank==0 && inf_local==numDirectionLocal-1 && iband_local==numBandLocal-1) {
                    cout<<"Finish calculating coefficient"<<endl;
                }

                Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double>> cg1;
                sol1=cg1.compute(stiffMatrix).solve(Re1);
                for (int i = 0; i < numCell; ++i)
                {
                    energyDensity[iband_local][inf_local][i]= sol1(i) ;
                }
                auto solve_end = chrono::high_resolution_clock::now();
                solver1_time += chrono::duration_cast<chrono::microseconds>(solve_end - solve_start);

                auto boundary_start = chrono::high_resolution_clock::now();
                _get_bound_ee_firstorder(iband_local,inf_local);
                auto boundary_end = chrono::high_resolution_clock::now();
                Boundary_time+=chrono::duration_cast<chrono::microseconds>(boundary_end - boundary_start);

                auto transfer_start = chrono::high_resolution_clock::now();
                _recover_temperature(iband_local,inf_local);
                _get_total_energy(iband_local,inf_local);
                _get_heat_flux(iband_local,inf_local);


                MPI_Allreduce(totalEnergyLocal, totalEnergy, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(temperatureLocal, temperature, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxXLocal, heatFluxXGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxYLocal, heatFluxYGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(heatFluxZLocal, heatFluxZGlobal, numCell, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                auto transfer_end = chrono::high_resolution_clock::now();
                trasfer1_time+=chrono::duration_cast<chrono::microseconds>(transfer_end - transfer_start);
            }
        }

        auto set_bound_start=chrono::high_resolution_clock::now();
        _set_bound_ee_1();
        auto set_bound_end=chrono::high_resolution_clock::now();
        set_bound_time+=chrono::duration_cast<chrono::microseconds>(set_bound_end - set_bound_start);


        auto total_iter_end = chrono::high_resolution_clock::now();
        total_iter_time += chrono::duration_cast<chrono::microseconds>(total_iter_end - total_iter_start);
        MPI_Barrier(MPI_COMM_WORLD);

        if (_get_magin_check_error(nt, error_temp_limit, error_flux_limit))
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
        }
        if (nt >= 5000 && nt != Num_Max_Iter && nt % 100 == 0)
        {
            if (worldRank == 0)
                cout << "More than 5,000 iterations have been performed. Please try \"Synthetic\" in IterativeScheme" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (errorIncreaseTime >= 10)
        {
            nt = Num_Max_Iter;
            if (worldRank == 0)
                _print_out();
            if (worldRank == 0)
                cout << "Error increases for 10 times, maybe the solution can not converge. Try to set Order=1" << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);


    _get_bound_temp();
    _get_bound_flux();


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
    _delete_matrix("BICGSTAB");



 }



