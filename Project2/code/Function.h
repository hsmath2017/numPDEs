/**
 * @file   Function.h
 * @author Shuang Hu <hsmath@ubuntu>
 * @date   Thu Jun  3 05:14:37 2021
 * 
 * @brief  Class Function.h:Mark the right-hand of the BVP
 * 
 * 
 */
#include<iostream>
#include<cmath>
double func(double x){
    double y=sin(x)*exp(sin(x))-pow(cos(x),2)*exp(sin(x));
    return y;
}


