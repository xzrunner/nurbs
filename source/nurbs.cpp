#include "../include/nurbs/nurbs.h"
#include "../external/aitn/bezier.h"

#include <memory>

namespace nurbs
{

void bezier(const std::vector<sm::vec2>& ctl_pts, std::vector<sm::vec2>& polyline)
{
	if (ctl_pts.empty() || polyline.size() < 2) {
		return;
	}

	auto dump_vec2 = [](const std::vector<sm::vec2>& array) -> std::vector<float>
	{
		std::vector<float> ret;
		ret.resize(array.size() * 2);
		for (int i = 0, n = array.size(); i < n; ++i) 
		{
			auto& p = array[i];
			ret[i * 2]     = p.x;
			ret[i * 2 + 1] = p.y;
		}
		return ret;
	};

	auto b = dump_vec2(ctl_pts);
	auto p = dump_vec2(polyline);
	nlib::bezier<float, 2>(ctl_pts.size(), b.data(), polyline.size(), p.data());

	for (int i = 0, n = polyline.size(); i < n; ++i) 
	{
		polyline[i].x = p[i * 2];
		polyline[i].y = p[i * 2 + 1];
	}
}

}