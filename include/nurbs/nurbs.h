#pragma once

#include <../sm/SM_Vector.h>

#include <vector>

namespace nurbs
{

void bezier(const std::vector<sm::vec2>& ctl_pts, 
	std::vector<sm::vec2>& polyline);

}