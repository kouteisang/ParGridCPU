

#include "Graph/MultilayerGraph.h"
#include "Core/FCTree.h"
#include "header.h"
#include "Util/MemoryUtils.h"
#include "Util/KLmdUtils.h"
#include "Util/CheckMLCDTime.h"
#include "Util/CoreIndexTime.h"

#include <iomanip>
// Naive
#include "Core/FCTreeDFS.h"

// serial
#include "Core/FCTreeBuilderRight.h"
#include "Core/FCTreeBuilderLeft.h"
#include "Core/FCNaive.h"
// new path level
#include "Core/FCPathLevelLeft.h"
#include "Core/FCPathLevelRight.h"

// path level
#include "Core/FCTreeBuilderPathParallelByk.h"
#include "Core/FCTreeBuilderPathParallelBylmd.h"

// core level
#include "CoreParallel/FCCoreTree.cpp"
#include "CoreParallel/FCTreeBuilderCoreParallel.h"
#include "CoreParallel/FCTreeBuilderCoreParallelByK.h"
#include "CoreParallel/FCSyncLeft.h"
#include "CoreParallel/FCSyncRight.h"

// Memory reduction
#include "Core/FCPathMem.h"
#include "CoreParallel/FCSyncMem.h"

// A new test
#include "CoreParallelNew/CoreParallelNew.h"

// A new algorithm
#include "CoreIndexParallel/CoreIndex.h"



uint* getLmd(){

}

int main(int argc, char* argv[]){
     

    string dataset = "example";
    string method = "serial";
    int order = 0;
    uint k = 0;
    uint lmd = 0;
    uint num_thread = 1;
    double beta = 1.0;

    for(int i = 1; i < argc; i ++){
        string arg = argv[i];
        if (arg == "-d" && i + 1 < argc) {
            dataset = argv[++i];
        }
        if(arg == "-m" && i + 1 < argc){
            method = argv[++i];
        }
        if(arg == "-o" && i + 1 < argc){
            order = std::atoi(argv[++i]);
        }
        if(arg == "-k" && i + 1 < argc){
            k = std::atoi(argv[++i]);
        }
        if(arg == "-lmd" && i + 1 < argc){
            lmd = std::atoi(argv[++i]);
        }
        if(arg == "-num_thread" && i+1 < argc){
            num_thread = std::atoi(argv[++i]);
        }
        if(arg == "-beta" && i+1 < argc){
            beta = std::atof(argv[++i]);
        }
    }

    if(dataset == ""){
        cout << "Use the default dataset example as the user did not specify the dataset" << endl;
    }
    if(method == ""){
        cout << "Use the FirmCore Tree Serial method as the user did not specify the method" << endl;
    }

    cout << "method = " << method  << " dataset = " << dataset << endl;

    // load the dataset
    MultilayerGraph mg;
    // relative path is relative to the execute file
    mg.LoadFromFile("../dataset/"+dataset+"/");
    mg.SetGraphOrder(order);
    uint *order_list = mg.GetOrder();
    for(int i = 0; i < mg.getLayerNumber(); i ++){
        cout << order_list[i] << " ";;
    }
    cout << endl;
    mg.PrintStatistics();

   
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber();

    if(method == "fcgrid"){
        uint lmdupper = mg.getLayerNumber();
        int kupper = 0;

        for(int l = 0; l < n_layers; l ++){
            int cnt = 0;
            int maxx = 0;
            for(int v = 0; v < n_vertex; v ++){
                int d = mg.GetGraph(l).GetAdjLst()[v][0];
                maxx = std::max(maxx, d);
            }
            kupper = std::max(kupper, maxx); 
        }
        if(lmdupper <= kupper){
            method = "lmdlist";
        }else{
            method = "klist";
        }
        cout << "kupper = " << kupper << endl;
        cout << "lmdupper = " << lmdupper << endl;
        cout << "method = " << method << endl;
    }

    ll_uint *id2vtx = new ll_uint[mg.GetN()];
    mg.LoadId2VtxMap(id2vtx);

    // delete
    if(method == "naive"){
        auto start_time = omp_get_wtime();
        FCTreeDFS::Execute(mg);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "DFS Elapsed time: " << elapsed_time << " seconds\n";
        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl;
 
    }

    // reviewer baseline
    if(method == "naivekl"){

        auto core_path = "/home/cheng/MlcDec/leaf/"+ dataset +"_cores.txt";
        auto klmd_path = "/home/cheng/fctree/klmd/"+ dataset +"_klmd.txt";

        vector<int> ks;
        vector<int> lmds;
        
        vector<vector<int>> cores;

        std::ifstream inputFile(klmd_path); 
        int kk, lmdd;
        while (inputFile >> kk >> lmdd) {
            ks.push_back(kk);
            lmds.push_back(lmdd);
        }

        auto start_time = omp_get_wtime();

        for(int p = 0; p < ks.size(); p ++){
            FCNaive::Execute(mg, ks[p], lmds[p]);
        }

        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "FC Naive Elapsed time: " << elapsed_time << " seconds\n";
        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl;
 
    }

    // Serial
    if(method == "lpe"){
        auto start_time = omp_get_wtime(); 
        FCTree tree(1, 1, mg.GetN());
        FCTreeBuilderLeft::Execute(mg, tree);
        // FCTreeBuilderLeft::Execute(mg);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "Serial Elapsed time: " << elapsed_time << " seconds\n";
        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl; 
    }

    if(method == "kpe"){
        auto start_time = omp_get_wtime();
        // This method simply traverse the tree
        // FCTreeBuilder::Execute(mg);

        // This method get the result and store the result in the tree
        FCTree tree(1, 1, mg.GetN());
        FCTreeBuilderRight::Execute(mg, tree);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "Serial Elapsed time: " << elapsed_time << " seconds\n";
        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl;

        // This part to get the output
        if(k != 0 && lmd != 0){
            k = uint(k);
            lmd = uint(lmd);
            coreNode* res_node = tree.getCoreByKAndLmdByLeft(tree.getNode(), k, lmd);
            cout << k << " " << lmd << endl;
            if(res_node == nullptr){
                cout << "No statisfy result" << endl;
            }else{
                tree.saveCoreToLocal(dataset, id2vtx, res_node, "serial");
            }
        }

    }



    // Core Parallel
    if(method == "CoreParallel"){
        omp_set_num_threads(num_thread);
        auto start_time = omp_get_wtime(); 
        FCCoreTree tree(1, 1, mg.GetN());
        // coreNodeP* node = tree.getNode();
        FCTreeBuilderCoreParallel::Execute(mg, tree);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "Path left Core Parallel Elapsed time: " << elapsed_time << " seconds\n";

        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl;
    }

    // Core Parallel with its own buffer
    if(method == "ParGrid"){
        omp_set_num_threads(num_thread);
        cout << "#thread = " << num_thread << endl; 
        auto start_time = omp_get_wtime(); 
        FCCoreTree tree(1, 1, mg.GetN());
        // coreNodeP* node = tree.getNode();
        FCSyncLeft::Execute(mg, tree);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "Path left Core Parallel Sync Elapsed time: " << elapsed_time << " seconds\n";

        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl; 

        string resfilename = dataset + "-" + method + "-" + std::to_string(num_thread);
        std::ofstream outfile(resfilename, std::ios::app); // 追加模式
        if (outfile.is_open()) {
            outfile << std::fixed << std::setprecision(6)
                    << "Elapsed time: " << elapsed_time << " seconds\n"
                    << "Memory usage: " << mem << " MB\n"
                    << "-----------------------------\n";
            outfile.close();
            std::cout << "Results saved to result.txt\n";
        } else {
            std::cerr << "Error: unable to open result.txt\n";
        }
        
    }

        // Core Parallel with its own buffer
    if(method == "CoreParallelSyncRight"){
        omp_set_num_threads(num_thread);
        auto start_time = omp_get_wtime(); 
        FCCoreTree tree(1, 1, mg.GetN());
        // coreNodeP* node = tree.getNode();
        FCSyncRight::Execute(mg, tree);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "Path Right Core Parallel Sync Elapsed time: " << elapsed_time << " seconds\n";

        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl; 
    }
    

    // Mix strategy with Core Parallel and Path Parallel
    if(method == "ParGridPlus"){
        omp_set_nested(1); // allow nest parallelization
        omp_set_max_active_levels(2); // maximum 2 layers
       
        omp_set_num_threads(num_thread);

        auto start_time = omp_get_wtime(); 
        FCCoreTree tree(1, 1, mg.GetN());
        // coreNodeP* node = tree.getNode();
        FCSyncLeft::ExecuteMix(mg, tree);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "Mix Core/Path Parallel Sync Elapsed time: " << elapsed_time << " seconds\n";

        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl; 

        string resfilename = dataset + "-2-" + method + "-" + std::to_string(num_thread);
        std::ofstream outfile(resfilename, std::ios::app); // 追加模式
        if (outfile.is_open()) {
            outfile << std::fixed << std::setprecision(6)
                    << "Elapsed time: " << elapsed_time << " seconds\n"
                    << "Memory usage: " << mem << " MB\n"
                    << "-----------------------------\n";
            outfile.close();
            std::cout << "Results saved to result.txt\n";
        } else {
            std::cerr << "Error: unable to open result.txt\n";
        }
    }


    // Mix strategy with Core Parallel and Path Parallel
    if(method == "ParGridPlusright"){
        omp_set_nested(1); // 允许嵌套并行
        omp_set_max_active_levels(2); // 最多允许两层并行
        
        omp_set_num_threads(num_thread);

        auto start_time = omp_get_wtime(); 
        FCCoreTree tree(1, 1, mg.GetN());
        // coreNodeP* node = tree.getNode();
        FCSyncRight::ExecuteMix(mg, tree);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "Mix Core/Path Parallel Sync Elapsed time: " << elapsed_time << " seconds\n";

        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl; 
    }

    // CoreParallel with memory reduction strategy
    if(method == "CoreMem"){
        omp_set_num_threads(num_thread);

        auto start_time = omp_get_wtime(); 
        FCCoreTree tree(1, 1, mg.GetN());
        // coreNodeP* node = tree.getNode();
        FCSyncMem::Execute(mg, tree);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "CoreParallel Sync with Memory reduction Elapsed time: " << elapsed_time << " seconds\n";

        long double mem = GetPeakRSSInMB();
        cout << "CoreParallel Sync with Memory reduction mem = " << mem << " MB" << endl; 
    }
    
    if(method == "PathMem"){
        omp_set_num_threads(num_thread);

        auto start_time = omp_get_wtime(); 
        FCTree tree(1, 1, mg.GetN());
        coreNode* node = tree.getNode();
        FCPathMem::Execute(mg, tree);
        auto end_time = omp_get_wtime(); 
        
        double elapsed_time = end_time - start_time;
        std::cout << "PathParallel with Memory reduction Elapsed time: " << elapsed_time << " seconds\n";

        long double mem = GetPeakRSSInMB();
        cout << "PathParallel Sync with Memory reduction mem = " << mem << " MB" << endl;  
    }



    // Baseline

    if(method == "CoreIndex"){
      
        auto start_time = omp_get_wtime();
        CoreIndex::Execute(mg, id2vtx);
        auto end_time = omp_get_wtime(); 
        double elapsed_time = end_time - start_time;
        std::cout << "CoreIndex Serial Elapsed time: " << elapsed_time << " seconds\n"; 
        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl;   
    }

    if(method == "ParCoreIndex"){
        omp_set_num_threads(num_thread);
        auto start_time = omp_get_wtime();
        CoreIndex::ExecuteParallel(mg, id2vtx);
        auto end_time = omp_get_wtime(); 
        double elapsed_time = end_time - start_time;
        std::cout << "CoreIndex Parallel Elapsed time: " << elapsed_time << " seconds\n"; 
        long double mem = GetPeakRSSInMB();
        cout << "mem = " << mem << " MB" << endl;

        string resfilename = dataset + "-" + method + "-" + std::to_string(num_thread);
        std::ofstream outfile(resfilename, std::ios::app); // 追加模式
        if (outfile.is_open()) {
            outfile << std::fixed << std::setprecision(6)
                    << "Elapsed time: " << elapsed_time << " seconds\n"
                    << "Memory usage: " << mem << " MB\n"
                    << "-----------------------------\n";
            outfile.close();
            std::cout << "Results saved to result.txt\n";
        } else {
            std::cerr << "Error: unable to open result.txt\n";
        }
    }

    if(method == "WdsList"){
        uint res_k = 0;
        uint res_lmd = 0;
        uint res_len = 0;
        float maximum_density = 0.0f; 
        omp_set_num_threads(num_thread);
        auto start_time = omp_get_wtime();
        FCTree tree(1, 1, mg.GetN());
        FCPathLevelLeft::Execute(mg, tree, true, beta);
        auto end_time = omp_get_wtime(); 


        #pragma omp parallel
        {
            #pragma omp single
            {
                tree.WeightDenestSubgraph(tree.getNode(), maximum_density, res_k, res_lmd, res_len);
            }
        }
        cout << "maximum_density = " << maximum_density << endl;
        cout << "res_k = " << res_k << endl;
        cout << "res_lmd = " << res_lmd << endl;
        cout << "res_len = " << res_len << endl;
         

        double elapsed_time = end_time - start_time;
        std::cout << "Wds Fctree Elapsed time: " << elapsed_time << " seconds\n";

        string resfilename = dataset + "-" + method + "-" + std::to_string(num_thread);
        std::ofstream outfile(resfilename, std::ios::app); // 追加模式
        if (outfile.is_open()) {
            outfile << std::fixed << std::setprecision(6)
                    << "Elapsed time: " << elapsed_time << " seconds\n"
                    << "-----------------------------\n";
            outfile.close();
            std::cout << "Results saved to result.txt\n";
        } else {
            std::cerr << "Error: unable to open result.txt\n";
        }

    }

    if(method == "WdsCoreIndex"){
        omp_set_num_threads(num_thread);
        auto start_time = omp_get_wtime();
        CoreIndex::WdsCoreIndex(mg, id2vtx, beta);
        auto end_time = omp_get_wtime(); 
        double elapsed_time = end_time - start_time;
        cout << "Wds CoreIndex Time = " << elapsed_time << endl;
    }

    if(method == "MLCDTime"){
        auto core_path = "/home/cheng/MlcDec/leaf/"+ dataset +"_cores.txt";
        auto klmd_path = "/home/cheng/fctree/klmd/"+ dataset +"_klmd.txt";

        vector<int> ks;
        vector<int> lmds;
        
        vector<vector<int>> cores;

        std::ifstream inputFile(klmd_path); 
        int kk, lmdd;
        while (inputFile >> kk >> lmdd) {
            ks.push_back(kk);
            lmds.push_back(lmdd);
        }

        string line;
        std::ifstream coreFile(core_path); 
        while (std::getline(coreFile, line)) {
            std::istringstream stream(line);  // 创建字符串流用于拆分
            std::vector<int> row;             // 存储一行的整数
            int number;

            // 将一行中的整数依次提取到row中
            while (stream >> number) {
                row.push_back(number);
            }

            // 将完整的一行（vector<int>）添加到data中
            cores.push_back(row);
        }

        auto start_time = omp_get_wtime();   
        checktime(ks, lmds, cores);
        auto end_time = omp_get_wtime();
        double elapsed_time = end_time - start_time;
        std::cout << "Elapsed time: " << elapsed_time << " seconds\n"; 

    }

    if(method == "CoreIndexTimeOnline"){

        auto klmd_path = "/home/cheng/fctree/klmd/"+ dataset +"_klmd.txt";

        vector<vector<int>> klmd(mg.getLayerNumber()+1);
        
        vector<vector<int>> cores;

        std::ifstream inputFile(klmd_path); 
        int kk, lmdd;
        while (inputFile >> kk >> lmdd) {
            klmd[lmdd].push_back(kk);
        }

        vector<vector<pair<int, int>>> results(mg.getLayerNumber()+1);
        int key, value;

        for(uint l = 1; l <= 1; l ++){
            std::ifstream coreFile("/home/cheng/fctree/Output/"+dataset+"-"+std::to_string(l)+".txt");

            while (coreFile >> key) {           // 读取第一个整数
                coreFile.ignore(1, ':');        // 忽略冒号
                coreFile >> std::ws >> value;    // 忽略空格并读取第二个整数
                results[l].push_back(std::make_pair(key, value)); // 存入pair并加入vector
            }

            coreFile.close();
        }
        auto start_time = omp_get_wtime(); 

        getResCore(klmd, results, mg.getLayerNumber());

        auto end_time = omp_get_wtime();  
        
        double retrieve_time = end_time - start_time;
        std::cout << "Elapsed time: " << retrieve_time << " seconds\n";
    }



    if(method == "CoreIndexTimeOffline"){

        auto klmd_path = "/home/cheng/fctree/klmd/"+ dataset +"_klmd.txt";

        vector<vector<int>> klmd(mg.getLayerNumber()+1);
        
        vector<vector<int>> cores;

        std::ifstream inputFile(klmd_path); 
        int kk, lmdd;
        while (inputFile >> kk >> lmdd) {
            klmd[lmdd].push_back(kk);
        }

        vector<vector<pair<int, int>>> results(mg.getLayerNumber()+1);
        int key, value;

        for(uint l = 1; l <= 1; l ++){
            std::ifstream coreFile("/home/cheng/fctree/Output/"+dataset+"-"+std::to_string(l)+".txt");

            while (coreFile >> key) {           // 读取第一个整数
                coreFile.ignore(1, ':');        // 忽略冒号
                coreFile >> std::ws >> value;    // 忽略空格并读取第二个整数
                results[l].push_back(std::make_pair(value, key)); // 存入pair并加入vector
            }

            coreFile.close();
        }
        auto start_time = omp_get_wtime(); 

        for (auto& vec : results) {
            std::sort(vec.begin(), vec.end(),
                    [](const std::pair<int,int>& a, const std::pair<int,int>& b){
                        return a.first < b.first;
                    });
        }
        
        for(uint l = 1; l <= 1; l ++){
            std::vector<int> ks = klmd[l];

            for(uint i = 0; i < ks.size(); i ++){
                uint k = ks[i];
                getResCoreOffline(klmd, results, mg.getLayerNumber(), k, l);
            }

        }

        auto end_time = omp_get_wtime();  
        
        double retrieve_time = end_time - start_time;
        std::cout << "Elapsed time: " << retrieve_time << " seconds\n";
    }




 // Path Parallel [[deprecated]]
    // if(method == "ParList"){
    //     omp_set_num_threads(num_thread);
    //     auto start_time = omp_get_wtime();
    //     FCTree tree(1, 1, mg.GetN());
    //     FCPathLevelLeft::Execute(mg, tree, false, beta);
    //     auto end_time = omp_get_wtime(); 
        
    //     double elapsed_time = end_time - start_time;
    //     std::cout << "PathLevelLeft Elapsed time: " << elapsed_time << " seconds\n";
    //     long double mem = GetPeakRSSInMB();
    //     cout << "mem = " << mem << " MB" << endl;
    //      // This part to get the output
    //     // std::vector<int> ks = getK("/home/cheng/fctree/dataset/s1/k.txt");
    //     // std::vector<int> lmds = getLmd("/home/cheng/fctree/dataset/s1/lmd.txt");

   
    //     // auto start_time_query = omp_get_wtime();
    //     // for(int i = 0; i < ks.size(); i ++){
    //     //     uint k = uint(ks[i]);
    //     //     uint lmd = uint(lmds[i]); 
    //     //     coreNode* res_node = tree.getCoreByKAndLmdByRight(tree.getNode(), k, lmd);
    //     //     // cout << res_node->k << " " << res_node->lmd << " " << res_node->length << endl;
    //     // }
    //     // auto end_time_query = omp_get_wtime();
    //     // auto total_time_query = end_time_query - start_time_query;
    //     // cout << "total time query = "<< total_time_query << endl;
    //     // cout << "average time query = "<< total_time_query/ks.size() << endl;
        
    //     if(k != 0 && lmd != 0){
    //         k = uint(k);
    //         lmd = uint(lmd);
    //         coreNode* res_node = tree.getCoreByKAndLmdByRight(tree.getNode(), k, lmd);
    //         if(res_node == nullptr){
    //             cout << "No statisfy result" << endl;
    //         }else{
    //             tree.saveCoreToLocal(dataset, id2vtx, res_node, "pathlmd");
    //         }
    //     }
    // }


//    // Path Parallel[[deprecated]]
//    if(method == "ParListlRight"){
//         omp_set_num_threads(num_thread);
//         auto start_time = omp_get_wtime();
//         FCTree tree(1, 1, mg.GetN());
//         FCPathLevelRight::Execute(mg, tree);
//         auto end_time = omp_get_wtime(); 
        
//         double elapsed_time = end_time - start_time;
//         std::cout << "PathLevelRight Elapsed time: " << elapsed_time << " seconds\n";
//         long double mem = GetPeakRSSInMB();
//         cout << "mem = " << mem << " MB" << endl;
//     }

}