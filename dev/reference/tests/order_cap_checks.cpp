#include "FunctionCall.h"
#include "WeightedMomentCore.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace wmc = weighted_moment_core;
static double reference(const std::vector<double>& values, int order,
                        const std::vector<double>& weights = std::vector<double>()) {
    long double mass=0, mean=0;
    for(size_t i=0;i<values.size();++i) {
        double w=weights.empty()?1:weights[i]; mass+=w; mean+=w*values[i];
    }
    mean/=mass;
    std::vector<long double> m(order+1,0), k(order+1,0); m[0]=1;
    for(int r=1;r<=order;++r) {
        for(size_t i=0;i<values.size();++i) {
            long double p=1; for(int j=0;j<r;++j)p*=values[i]-mean;
            m[r]+=(weights.empty()?1:weights[i])*p/mass;
        }
        k[r]=m[r]; long double choose=1;
        for(int j=1;j<r;++j) { k[r]-=choose*k[j]*m[r-j]; choose=choose*(r-j)/j; }
    }
    return order==1?(double)mean:(double)k[order];
}
static void close(double a,double b) {
    if(!std::isfinite(a) || std::abs(a-b)>1e-12+2e-8*std::abs(b))
        throw std::runtime_error("High-order numerical mismatch: "+std::to_string(a)+" vs "+std::to_string(b));
}
static double avl(const std::vector<double>& values,const std::vector<int>& powers) {
    int n=(int)powers.size(); std::vector<int> index(n+1,0);
    for(int j=0;j<n;++j){index[j+1]=powers[j];index[0]+=powers[j];}
    double mean=0;for(double v:values)mean+=v/values.size();
    std::vector<std::vector<double>> rows(values.size(),std::vector<double>(n+1,0));
    std::vector<Template> tuples(values.size());
    for(size_t i=0;i<values.size();++i) {
        for(int j=1;j<=n;++j)rows[i][j]=values[i]-mean;
        tuples[i].Value=rows[i].data(); tuples[i].index=index.data();
    }
    ClearOneEvaluationWorkspaceForDLL(n);
    G_Template=tuples.data(); G_TemplateNum=(int)tuples.size();
    try {
        SetFactorSwitchesForRun(1,1,1,1);SetLowOrderCumulantMomentCacheMode(1);
        if(PrecomputeMomentsForIndex(index.data(),n)==myEmpty || InitAVLWorkspaceForDLL(n)!=0)
            throw std::runtime_error("Failed to initialize AVL test.");
        double result=BetterRecCumulant3(index.data(),n);
        ClearOneEvaluationWorkspaceForDLL(n);G_Template=NULL;G_TemplateNum=0;
        return result;
    } catch(...) {
        ClearOneEvaluationWorkspaceForDLL(n);G_Template=NULL;G_TemplateNum=0;throw;
    }
}
template<class F> static void rejects(F f,const char* text) {
    try {f();} catch(const std::exception& e) {
        if(std::string(e.what()).find(text)!=std::string::npos)return;
        throw;
    }
    throw std::runtime_error("Expected explicit capacity failure.");
}
int main() {
    try {
        std::vector<double> v{1.1,1.21,1.3,1.44,1.53,1.81};
        const std::vector<std::vector<int>> indices{{1,4,4},{1,4,5},{1,5,5},{1,5,6},{1,8,7}};
        for(const auto& powers:indices) {
            int order=0;for(int p:powers)order+=p;
            close(avl(v,powers),reference(v,order));
            std::cout<<"AVL template order "<<order<<" passed\n";
        }
        rejects([&]{avl(v,{10});},"existing exponent encoding");
        int index[]={10,1,4,5};
        rejects([&]{BellCumulantFromMultiplicityIndex(index,3);},"workspace capacity");
        std::vector<double> weights{1,2,3,4,2,1};
        for(int order:{10,12,16})
            close(wmc::weighted_univariate_cumulant(v.data(),NULL,(int)v.size(),order,true,true),reference(v,order));
        wmc::WeightedMomentOptions o; o.mode=wmc::WeightMode::External; o.zero_first_order_moment=true;
        std::vector<std::vector<double>> rows(v.size(),std::vector<double>(3));
        std::vector<wmc::WeightedSample> samples(v.size());
        for(size_t i=0;i<v.size();++i){rows[i].assign(3,v[i]);samples[i].values=rows[i].data();samples[i].value_count=3;samples[i].weight=weights[i];}
        int powers[]={1,4,5};
        close(wmc::weighted_cumulant_from_multiplicities(samples.data(),(int)v.size(),powers,3,o),reference(v,10,weights));
        std::cout<<"Weighted Bell template order 10 and window recurrence orders 10/12/16 passed\n";
        powers[2]=27;
        rejects([&]{wmc::weighted_moment(samples.data(),(int)v.size(),powers,3,o);},"partition-mask capacity");
        std::cout<<"Unchanged core capacity failures are explicit; no truncated results\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
