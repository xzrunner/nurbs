#pragma once

#include <../sm/SM_Vector.h>

#include <vector>

namespace nurbs
{

void bezier(const std::vector<sm::vec2>& ctl_pts, 
	std::vector<sm::vec2>& polyline);
void bspline(const std::vector<sm::vec2>& ctl_pts,
	int order, std::vector<sm::vec2>& polyline);
void rbspline(const std::vector<sm::vec2>& ctl_pts,
	int order, std::vector<sm::vec2>& polyline);
void bezsurf(const sm::vec3* ctl_pts, int npts, int mpts,
	int p1, int p2, std::vector<sm::vec3>& surface);
void bspsurf(const sm::vec3* ctl_pts, int order_u, int order_v,
	int npts, int mpts, int p1, int p2, std::vector<sm::vec3>& surface);
void rbspsurf(const sm::vec3* ctl_pts, int order_u, int order_v,
	int npts, int mpts, int p1, int p2, std::vector<sm::vec3>& surface);

}