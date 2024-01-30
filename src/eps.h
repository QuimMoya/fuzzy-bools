#pragma once

static double EPS_SCALE = 1;

static constexpr double E_NONZERO = 1e-20; // ???
static constexpr double E_MINISCULE = 1e-12; // what?
static constexpr double E_TINY = 1e-9;
static constexpr double E_SMALL = 1e-6;
static constexpr double E_BIG = 1e-4;
static constexpr double E_BIG2 = 1e-3;
static constexpr double E_BIG3 = 1e-2;

static double RET_VAL(double val, bool scale)
{
    if (scale) 
    {
        return val / EPS_SCALE;
    }
    else
    {
        return val;
    }
}

static double EPS_NONZERO(bool scale = false){
    return RET_VAL(E_NONZERO, scale);
}

static double EPS_MINISCULE(bool scale = false){
    return RET_VAL(E_MINISCULE, scale);
}

static double EPS_TINY(bool scale = false){
    return RET_VAL(E_TINY, scale);
}

static double EPS_SMALL(bool scale = false){
    return RET_VAL(E_SMALL, scale);
}

static double EPS_BIG(bool scale = false){
    return RET_VAL(E_BIG, scale);
}

static double EPS_BIG2(bool scale = false){
    return RET_VAL(E_BIG2, scale);
}

static double EPS_BIG3(bool scale = false){
    return RET_VAL(E_BIG3, scale);
}