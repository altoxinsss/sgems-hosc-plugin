#ifndef HOSC_RUN_PERFORMANCE_H
#define HOSC_RUN_PERFORMANCE_H

// Per-execute host metrics. No numerical-core state is read or changed.
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef PSAPI_VERSION
#define PSAPI_VERSION 1
#endif
#include <windows.h>
#include <psapi.h>
#include <process.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <map>
#include <algorithm>
#pragma comment(lib,"psapi.lib")

namespace hosc_sgems {
class RunPerformance {
public:
    enum { sample_interval_ms = 5 };
    RunPerformance():event_(0),thread_(0),stopped_(false),timer_valid_(false),
        baseline_(0),peak_(0),private_peak_(0),lifetime_peak_(0),samples_(0),failures_(0),sampler_started_(false) {
        InitializeCriticalSection(&mutex_);
        ZeroMemory(&frequency_,sizeof(frequency_));ZeroMemory(&start_,sizeof(start_));ZeroMemory(&end_,sizeof(end_));
        GetSystemTime(&started_);
        timer_valid_=QueryPerformanceFrequency(&frequency_)!=0 && frequency_.QuadPart>0 && QueryPerformanceCounter(&start_)!=0;
        end_=start_;
        sample(); baseline_=peak_;
        event_=CreateEventA(0,TRUE,FALSE,0);
        if(event_)thread_=reinterpret_cast<HANDLE>(_beginthreadex(0,0,&RunPerformance::worker,this,0,0));
        sampler_started_=thread_!=0;
    }
    ~RunPerformance(){stop();if(thread_)CloseHandle(thread_);if(event_)CloseHandle(event_);DeleteCriticalSection(&mutex_);}

    void checkpoint(){sample();}
    void stop(){
        if(stopped_)return;
        // Final synchronous reading covers even calculations shorter than 5 ms.
        sample();
        if(timer_valid_ && !QueryPerformanceCounter(&end_))timer_valid_=false;
        GetSystemTime(&finished_);
        if(event_)SetEvent(event_);
        if(thread_)WaitForSingleObject(thread_,INFINITE);
        stopped_=true;
    }
    double elapsed_seconds()const{return timer_valid_?double(end_.QuadPart-start_.QuadPart)/double(frequency_.QuadPart):0;}
    double baseline_mib()const{return double(baseline_)/1048576.0;}
    double peak_mib()const{return double(peak_)/1048576.0;}
    double lifetime_peak_mib()const{return double(lifetime_peak_)/1048576.0;}
    unsigned long samples()const{return samples_;}
    bool complete()const{return stopped_ && timer_valid_ && sampler_started_ && samples_>0 && failures_==0;}

    bool write(const std::string& directory,const std::map<std::string,std::string>& metadata,
               const std::string& status,const std::string& error) {
        stop();
        std::map<std::string,std::string> strings=metadata;
        strings["status"]=status;strings["error"]=error;
        strings["measurement_status"]=complete()?"ok":"incomplete";
        strings["started_utc"]=stamp(started_);strings["finished_utc"]=stamp(finished_);
        strings["run_directory"]=directory;
        strings["elapsed_scope"]="plugin execute: run-directory creation, input export, numerical calculation, diagnostics, result import and project update; excludes application startup, parameter loading and performance-record writing";
        strings["memory_scope"]="current SGeMS process working set observed during this execute call; includes host baseline, inputs and retained results";
        strings["memory_method"]="GetProcessMemoryInfo WorkingSetSize: 5-ms requested sampling interval plus phase checkpoints; sampled peak, not an exact instantaneous maximum";
        strings["repeat_role"]="unassigned";
        strings["application_executable"]=executable();
        std::map<std::string,double> numbers;
        numbers["schema_version"]=1;
        numbers["elapsed_seconds"]=elapsed_seconds();
        numbers["baseline_working_set_mib"]=baseline_mib();
        numbers["peak_working_set_mib"]=peak_mib();
        numbers["peak_working_set_increase_mib"]=double(peak_>=baseline_?peak_-baseline_:0)/1048576.0;
        numbers["peak_private_usage_mib"]=double(private_peak_)/1048576.0;
        // Kept only as an explicitly named diagnostic, never the per-run peak.
        numbers["process_lifetime_peak_working_set_mib"]=lifetime_peak_mib();
        numbers["memory_sample_count"]=samples_;
        numbers["memory_read_failures"]=failures_;
        numbers["memory_sampling_interval_ms"]=sample_interval_ms;
        numbers["memory_sampler_started"]=sampler_started_?1:0;
        numbers["timer_valid"]=timer_valid_?1:0;
        numbers["process_id"]=GetCurrentProcessId();
        numbers["process_bits"]=sizeof(void*)*8;
        std::ofstream json((directory+"\\performance.json").c_str());
        if(!json)return false;
        json<<std::setprecision(17)<<"{\n";bool first=true;
        for(std::map<std::string,std::string>::const_iterator i=strings.begin();i!=strings.end();++i){
            if(!first)json<<",\n";first=false;json<<"  \""<<escaped(i->first)<<"\": \""<<escaped(i->second)<<"\"";
        }
        for(std::map<std::string,double>::const_iterator i=numbers.begin();i!=numbers.end();++i){
            if(!first)json<<",\n";first=false;json<<"  \""<<i->first<<"\": ";
            if((i->first=="elapsed_seconds"&&!timer_valid_) ||
               ((i->first=="peak_working_set_mib"||i->first=="baseline_working_set_mib")&&samples_==0))json<<"null";
            else json<<i->second;
        }
        json<<"\n}\n";json.close();if(!json)return false;
        std::ofstream csv((directory+"\\performance.csv").c_str(),std::ios::out|std::ios::binary);if(!csv)return false;
        csv<<std::setprecision(17);first=true;
        for(std::map<std::string,std::string>::const_iterator i=strings.begin();i!=strings.end();++i){if(!first)csv<<',';first=false;csv<<quoted(i->first);}
        for(std::map<std::string,double>::const_iterator i=numbers.begin();i!=numbers.end();++i){if(!first)csv<<',';first=false;csv<<quoted(i->first);}
        csv<<'\n';first=true;
        for(std::map<std::string,std::string>::const_iterator i=strings.begin();i!=strings.end();++i){if(!first)csv<<',';first=false;csv<<quoted(i->second);}
        for(std::map<std::string,double>::const_iterator i=numbers.begin();i!=numbers.end();++i){
            if(!first)csv<<',';first=false;
            if(i->first=="elapsed_seconds"&&!timer_valid_)continue;
            if((i->first=="peak_working_set_mib"||i->first=="baseline_working_set_mib")&&samples_==0)continue;
            csv<<i->second;
        }
        csv<<'\n';csv.close();return bool(csv);
    }
private:
    RunPerformance(const RunPerformance&);
    RunPerformance& operator=(const RunPerformance&);
    HANDLE event_,thread_;CRITICAL_SECTION mutex_;bool stopped_,timer_valid_;
    LARGE_INTEGER frequency_,start_,end_;SYSTEMTIME started_,finished_;
    SIZE_T baseline_,peak_,private_peak_,lifetime_peak_;
    unsigned long samples_,failures_;bool sampler_started_;

    void sample(){
        EnterCriticalSection(&mutex_);
        PROCESS_MEMORY_COUNTERS_EX memory;ZeroMemory(&memory,sizeof(memory));memory.cb=sizeof(memory);
        if(GetProcessMemoryInfo(GetCurrentProcess(),reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&memory),sizeof(memory))){
            peak_=std::max(peak_,memory.WorkingSetSize);
            private_peak_=std::max(private_peak_,memory.PrivateUsage);
            lifetime_peak_=std::max(lifetime_peak_,memory.PeakWorkingSetSize);++samples_;
        }else ++failures_;
        LeaveCriticalSection(&mutex_);
    }
    static unsigned __stdcall worker(void* context){
        RunPerformance* self=static_cast<RunPerformance*>(context);
        while(WaitForSingleObject(self->event_,sample_interval_ms)==WAIT_TIMEOUT)self->sample();
        return 0;
    }
    static std::string stamp(const SYSTEMTIME& t){
        std::ostringstream s;s<<std::setfill('0')<<std::setw(4)<<t.wYear<<'-'<<std::setw(2)<<t.wMonth<<'-'<<std::setw(2)<<t.wDay<<'T'
            <<std::setw(2)<<t.wHour<<':'<<std::setw(2)<<t.wMinute<<':'<<std::setw(2)<<t.wSecond<<'.'<<std::setw(3)<<t.wMilliseconds<<'Z';return s.str();
    }
    static std::string executable(){char path[32768];DWORD n=GetModuleFileNameA(0,path,sizeof(path));return n&&n<sizeof(path)?std::string(path,n):"unavailable";}
    static std::string escaped(const std::string& s){
        std::ostringstream out;
        for(size_t i=0;i<s.size();++i){unsigned char c=static_cast<unsigned char>(s[i]);
            if(c=='"')out<<"\\\"";else if(c=='\\')out<<"\\\\";
            else if(c<32)out<<"\\u00"<<std::hex<<std::setfill('0')<<std::setw(2)<<int(c)<<std::dec;
            else out<<char(c);
        }return out.str();
    }
    static std::string quoted(const std::string& s){std::string out="\"";for(size_t i=0;i<s.size();++i){if(s[i]=='"')out+='"';out+=s[i];}return out+'"';}
};
}
#endif
