//
// Created by jsmith on 9/30/15.
//

#ifndef SIMULATOR_VECTOR_H
#define SIMULATOR_VECTOR_H

class point{
public:
    point(){
        x = 0;
        y = 0;
    }
    int x;
    int y;
};

class vec2{
public:
    vec2(){
        x = y = -1;
    }
    vec2(float _x, float _y){
        x = _x;
        y = _y;
    }
    vec2 operator=(const vec2& vector){
        if(this == &vector)
            return *this;
        x = vector.x;
        y = vector.y;
        return *this;
    }
    float x,y;
};
class vec3{
public:
    vec3(){
        x = y = z = -1;
    }
    vec3(float _x, float _y, float _z){
        x = _x;
        y = _y;
        z = _z;
    }
    vec3 operator=(const vec3& vector){
        if(this == &vector)
            return *this;
        x = vector.x;
        y = vector.y;
        z = vector.z;
        return *this;
    }
    float x,y,z;
};
class vec4{
public:
    vec4(){
        x = y = z = w = -1;
    }
    vec4(float _x, float _y, float _z, float _w){
        x = _x;
        y = _y;
        z = _z;
        w = _w;
    }
    vec4 operator=(const vec4& vector){
        if(this == &vector)
            return *this;
        x = vector.x;
        y = vector.y;
        z = vector.z;
        w = vector.w;
        return *this;
    }
    float x,y,z,w;
};

#endif //SIMULATOR_VECTOR_H
