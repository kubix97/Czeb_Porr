#pragma once
class Vector
{
    int     _iR;
    float*  _pv;

public:
    
    Vector() : _iR(0), _pv(nullptr)
    {}
    Vector(int length);
    Vector(const Vector& v);

    ~Vector();

    void                PrintVectorToShell();
    void                Generate(float min, float max);
    void                GenWithFixedVal(int len, float *vals);

    Vector              operator - (const Vector& v);
    Vector&             SubtractionI(Vector& v);

    Vector              Add(Vector& v);
    Vector&             AddI(Vector& v);

    Vector              MultiplyByVal(float value);
    Vector&             MultiplyByValI(float value);

    float               CalcDistance(Vector& v);

    int                 GetLen() const              { return _iR;       }
    float*              GetPtr() const              { return _pv;       }
    float               operator [] (int i) const   { return _pv[i];    }
    float&              operator [] (int i)         { return _pv[i];    }

    Vector&             operator = (const Vector& v);
};

