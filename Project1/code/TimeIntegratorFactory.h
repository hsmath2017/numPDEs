/**
 * @file   TimeIntegratorFactory.h
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Sat Apr 24 04:44:13 2021
 * 
 * @brief  Achieve the Factory Pattern
 * 
 * 
 */
#ifndef _HSMATH_FACTORY
#define _HSMATH_FACTORY
#include "AdamBashforce.h"
#include "BDF.h"
#include "AdamMoulton.h"
class TimeIntegratorFactory{
public:
    typedef TimeIntegrator* (*TimeIntegratorCreator)();
    bool RegisterTimeIntegrator(const string &TimeIntegratorID, TimeIntegratorCreator createTI);
    bool UnregisterTimeIntegrator(const string &TimeIntegratorID);
    TimeIntegrator* CreateTimeIntegrator(const string &TimeIntegratorID) const;
    static TimeIntegratorFactory* Instance(){
	if(!pInstance_)
	    pInstance_ = new TimeIntegratorFactory;
	return pInstance_;
    }
private:
    typedef map<string,TimeIntegratorCreator> CallbackMap;
    static TimeIntegratorFactory* pInstance_;
    TimeIntegratorFactory(){};
    TimeIntegratorFactory(const TimeIntegratorFactory&);
    CallbackMap callbacks_;
};
TimeIntegratorFactory* TimeIntegratorFactory::pInstance_ = 0;
namespace
{
    TimeIntegrator* CreateAM()
    {
	return new AdamMoulton;
    }
    TimeIntegrator* CreateAB()
    {
	return new AdamBashforce;
    }
    TimeIntegrator* CreateBDF()
    {
	return new BDF;
    }
    TimeIntegrator* CreateRK()
    {
	return new RungeKutta;
    }
    const bool _AM_registered = TimeIntegratorFactory::Instance()->RegisterTimeIntegrator("AdamMoulton",CreateAM);
    const bool _AB_registered = TimeIntegratorFactory::Instance()->RegisterTimeIntegrator("AdamBashforce",CreateAB);
    const bool _BDF_registered = TimeIntegratorFactory::Instance()->RegisterTimeIntegrator("BDF",CreateBDF);
    const bool _RK_registered = TimeIntegratorFactory::Instance()->RegisterTimeIntegrator("RungeKutta",CreateRK);
}
bool TimeIntegratorFactory::RegisterTimeIntegrator(const string &TimeIntegratorID,TimeIntegratorCreator createTI)
{
    return callbacks_.insert(
	CallbackMap::value_type(TimeIntegratorID,createTI)).second;
}
bool TimeIntegratorFactory::UnregisterTimeIntegrator(const string &TimeIntegratorID)
{
    return callbacks_.erase(TimeIntegratorID) ==1;
}
TimeIntegrator* TimeIntegratorFactory::CreateTimeIntegrator(const string &TimeIntegratorID) const
{
    auto i = callbacks_.find(TimeIntegratorID);
    if(i == callbacks_.cend()){//Not Found
	throw runtime_error("Unknown TimeIntegrator ID");
    }
    return (i->second)();
}
#else
//NOTHING
#endif
