#include "../src/shared/RunPerformance.h"
#include <stdexcept>
#include <iostream>

static void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(int argc,char** argv){
    try{
        require(argc==2,"output directory argument required");std::string out=argv[1];
        std::map<std::string,std::string> metadata;metadata["case_name"]="unit \"quote\", newline\ncheck";
        double first_peak=0;
        {
            hosc_sgems::RunPerformance metric;
            const SIZE_T bytes=64*1024*1024;
            volatile char* memory=static_cast<volatile char*>(VirtualAlloc(0,bytes,MEM_RESERVE|MEM_COMMIT,PAGE_READWRITE));
            require(memory!=0,"64 MiB allocation failed");
            for(SIZE_T i=0;i<bytes;i+=4096)memory[i]=1;
            Sleep(85);metric.checkpoint();
            VirtualFree(const_cast<char*>(memory),0,MEM_RELEASE);metric.stop();
            require(metric.complete(),"measurement incomplete");
            require(metric.elapsed_seconds()>=.075,"QPC elapsed time too small");
            require(metric.peak_mib()>=metric.baseline_mib()+48,"allocation not seen in working-set peak");
            require(metric.samples()>=3,"insufficient samples");first_peak=metric.peak_mib();
            CreateDirectoryA((out+"\\first").c_str(),0);
            require(metric.write(out+"\\first",metadata,"success",""),"record write failed");
        }
        Sleep(30);
        {
            hosc_sgems::RunPerformance metric;Sleep(45);metric.stop();
            require(metric.complete(),"second measurement incomplete");
            require(metric.peak_mib()<first_peak-32,"per-run peak incorrectly retained previous allocation");
            require(metric.lifetime_peak_mib()>metric.peak_mib()+32,"lifetime peak test did not distinguish histories");
            CreateDirectoryA((out+"\\second").c_str(),0);
            require(metric.write(out+"\\second",metadata,"failed","deliberate test error"),"failure record write failed");
        }
        DWORD handles_before=0,handles_after=0;
        require(GetProcessHandleCount(GetCurrentProcess(),&handles_before)!=0,"handle count failed");
        for(int i=0;i<20;++i){hosc_sgems::RunPerformance metric;metric.checkpoint();metric.stop();require(metric.complete(),"short run measurement failed");}
        require(GetProcessHandleCount(GetCurrentProcess(),&handles_after)!=0,"handle count failed");
        require(handles_after==handles_before,"monitor leaked process handles");
        std::cout<<"PASS timer, sampled allocation peak, reset versus lifetime peak, success/failure records, short runs and 20-run handle cleanup"<<std::endl;
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
