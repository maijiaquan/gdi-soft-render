#ifndef VERTEX_H_
#define VERTEX_H_

#include "Vector3.h"
#include "Math.h"

class Vertex
{
public:
	Vertex(const Vector3 &position,
		const AColor &color,
		float u,
		float v):
	position_(position),
		color_(color),
		u_(u),
		v_(v){}

	bool isxy_same(const Vertex &v)const{
		return position_.x_ == v.position_.x_&&position_.y_ == v.position_.y_;
	}
	Vertex interp(const Vertex &v,float factor);
	void interp_u_v(const Vertex &ver,float factor,float &u,float &v);
	Vector2 to_vector2(){ return Vector2(position_.x_,position_.y_,color_);}
	void set_normal(const Vector3 &normal){normal_ = normal;}

	Vector3 position_;
	Vector3 normal_;
	AColor color_;
	float u_,v_;
	Vector3 light_;
};


#endif