#include "../include/nurbs/nurbs.h"
#include "../external/aitn/bezier.h"
#include "../external/aitn/bspline.h"
#include "../external/aitn/rbspline.h"

#include <memory>

namespace
{

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

auto dump_vec3 = [](const std::vector<sm::vec2>& array) -> std::vector<float>
{
	std::vector<float> ret;
	ret.resize(array.size() * 3);
	for (int i = 0, n = array.size(); i < n; ++i) 
	{
		auto& p = array[i];
		ret[i * 3]     = p.x;
		ret[i * 3 + 1] = p.y;
		ret[i * 3 + 2] = 0;
	}
	return ret;
};

}

namespace nurbs
{

void bezier(const std::vector<sm::vec2>& ctl_pts, std::vector<sm::vec2>& polyline)
{
	if (ctl_pts.empty() || polyline.size() < 2) {
		return;
	}

	auto b = dump_vec2(ctl_pts);
	auto p = dump_vec2(polyline);
	aitn::bezier<float, 2>(ctl_pts.size(), b.data(), polyline.size(), p.data());

	for (int i = 0, n = polyline.size(); i < n; ++i) 
	{
		polyline[i].x = p[i * 2];
		polyline[i].y = p[i * 2 + 1];
	}
}

void bspline(const std::vector<sm::vec2>& ctl_pts, int order,
	         std::vector<sm::vec2>& polyline)
{
	if (ctl_pts.empty() || polyline.size() < 2) {
		return;
	}

	auto b = dump_vec3(ctl_pts);
	auto p = dump_vec3(polyline);
	aitn::bspline(ctl_pts.size(), order, polyline.size(), b.data(), p.data());
//	aitn::bsplineu(ctl_pts.size(), order, polyline.size(), b.data(), p.data());

	for (int i = 0, n = polyline.size(); i < n; ++i)
	{
		polyline[i].x = p[i * 3];
		polyline[i].y = p[i * 3 + 1];
	}
}

void rbspline(const std::vector<sm::vec2>& ctl_pts, int order,
	          std::vector<sm::vec2>& polyline)
{
	if (ctl_pts.empty() || polyline.size() < 2) {
		return;
	}

	std::vector<float> h(ctl_pts.size(), 1.0f);

	auto b = dump_vec3(ctl_pts);
	auto p = dump_vec3(polyline);
	aitn::rbspline(ctl_pts.size(), order, polyline.size(), b.data(), h.data(), p.data());
//	aitn::rbsplineu(ctl_pts.size(), order, polyline.size(), b.data(), h.data(), p.data());

	for (int i = 0, n = polyline.size(); i < n; ++i)
	{
		polyline[i].x = p[i * 3];
		polyline[i].y = p[i * 3 + 1];
	}
}

}