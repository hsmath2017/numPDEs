#ifndef _MG_Factory
#define _MG_Factory
#include "MGSolver.h"
#include <map>
template<int Dim>
class MGFactory{
public:
    using CreateMultigridSolverCallback = std::unique_ptr<MGSolver<Dim>>(*)(RectDomain<Dim>,const char*);
private:
    using CallbackMap = std::map<std::string,CreateMultigridSolverCallback>;
public:
    static MGFactory& createFactory(){
        static MGFactory object;
        return object;
    }

    bool registerMultigridSolver(std::string multigridId, CreateMultigridSolverCallback createFn){
        return callbacks__.insert(typename CallbackMap::value_type(multigridId,createFn)).second;
    }

    bool unregisterMultigridSolver(std::string multigirdId){
        return callbacks__.erase(multigirdId)==1;
    }

    template<class ...TS>
    std::unique_ptr<MGSolver<Dim>> CreateMultigridSolver(std::string multigridId, TS&& ...args){
        auto it=callbacks__.find(multigridId);
        if(it==callbacks__.end()){
            throw std::runtime_error("Unknown MultigridSolver ID. ");
        }
        return (it->second)(std::forward<TS>(args)...);
    }
    
private:
    MGFactory()=default;
    MGFactory(const MGFactory&)=default;
    MGFactory& operator=(const MGFactory&)=default;
    ~MGFactory()=default;
    CallbackMap callbacks__;
};
#else
//nothing
#endif