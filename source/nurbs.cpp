#include "../include/nurbs/nurbs.h"
#include "../external/aitn/bezier.h"
#include "../external/aitn/bspline.h"
#include "../external/aitn/rbspline.h"
#include "../external/aitn/bezsurf.h"
#include "../external/aitn/bspsurf.h"
#include "../external/aitn/rbspsurf.h"

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

void bezsurf(const sm::vec3* ctl_pts, int npts, int mpts,
	         int p1, int p2, std::vector<sm::vec3>& surface)
{
	float b[61];
	float q[301];

	memset(b, 0, sizeof(b));
	memset(q, 0, sizeof(q));

	b[1] = -15.;
	b[2] = 0.;
	b[3] = 15.;
	b[4] = -15.;
	b[5] = 5.;
	b[6] = 5.;
	b[7] = -15.;
	b[8] = 5.;
	b[9] = -5.;
	b[10] = -15.;
	b[11] = 0.;
	b[12] = -15.;

	b[13] = -5.;
	b[14] = 5.;
	b[15] = 15.;
	b[16] = -5.;
	b[17] = 5.;
	b[18] = 5.;
	b[19] = -5.;
	b[20] = 5.;
	b[21] = -5.;
	b[22] = -5.;
	b[23] = 5.;
	b[24] = -15.;

	b[25] = 5.;
	b[26] = 5.;
	b[27] = 15.;
	b[28] = 5.;
	b[29] = 5.;
	b[30] = 5.;
	b[31] = 5.;
	b[32] = 5.;
	b[33] = -5.;
	b[34] = 5.;
	b[35] = 5.;
	b[36] = -15.;

	b[37] = 15.;
	b[38] = 0.;
	b[39] = 15.;
	b[40] = 15.;
	b[41] = 5.;
	b[42] = 5.;
	b[43] = 15.;
	b[44] = 5.;
	b[45] = -5.;
	b[46] = 15.;
	b[47] = 0.;
	b[48] = -15.;

	for (int i = 0; i < 48; ++i) {
		b[i] = b[i + 1] * 0.01f;;
	}

	int n = npts - 1;
	int m = mpts - 1;
	//surface.resize(n * p1 * m * p2);
	//aitn::bezsurf(const_cast<float*>(ctl_pts[0].xyz), n, m, 
	//	p1, p2, surface[0].xyz);

	aitn::bezsurf(b, n, m, p1, p2, q);

	surface.resize(p1 * p2);

	int ptr = 0;
	for (int i = 0, n = p1 * p2; i < n; ++i) 
	{
		surface[i].x = q[ptr++];
		surface[i].y = q[ptr++];
		surface[i].z = q[ptr++];
	}
}

void bspsurf(const sm::vec3* ctl_pts, int order_u, int order_v,
	         int npts, int mpts, int p1, int p2, std::vector<sm::vec3>& surface)
{
	float b[61];
	float q[301];

	memset(b, 0, sizeof(b));
	memset(q, 0, sizeof(q));

	b[1] = -15.;
	b[2] = 0.;
	b[3] = 15.;
	b[4] = -15.;
	b[5] = 5.;
	b[6] = 5.;
	b[7] = -15.;
	b[8] = 5.;
	b[9] = -5.;
	b[10] = -15.;
	b[11] = 0.;
	b[12] = -15.;

	b[13] = -5.;
	b[14] = 5.;
	b[15] = 15.;
	b[16] = -5.;
	b[17] = 5.;
	b[18] = 5.;
	b[19] = -5.;
	b[20] = 5.;
	b[21] = -5.;
	b[22] = -5.;
	b[23] = 5.;
	b[24] = -15.;

	b[25] = 5.;
	b[26] = 5.;
	b[27] = 15.;
	b[28] = 5.;
	b[29] = 5.;
	b[30] = 5.;
	b[31] = 5.;
	b[32] = 5.;
	b[33] = -5.;
	b[34] = 5.;
	b[35] = 5.;
	b[36] = -15.;

	b[37] = 15.;
	b[38] = 0.;
	b[39] = 15.;
	b[40] = 15.;
	b[41] = 5.;
	b[42] = 5.;
	b[43] = 15.;
	b[44] = 5.;
	b[45] = -5.;
	b[46] = 15.;
	b[47] = 0.;
	b[48] = -15.;

	for (int i = 0; i < 48; ++i) {
		b[i] = b[i + 1] * 0.01f;;
	}

	int n = npts - 1;
	int m = mpts - 1;
	//surface.resize(n * p1 * m * p2);
	//aitn::bezsurf(const_cast<float*>(ctl_pts[0].xyz), n, m, 
	//	p1, p2, surface[0].xyz);

	aitn::bsplsurf(b, order_u, order_v, n, m, p1, p2, q);
	//aitn::bspsurfu(b, order_u, order_v, n, m, p1, p2, q);

	surface.resize(p1 * p2);

	int ptr = 0;
	for (int i = 0, n = p1 * p2; i < n; ++i) 
	{
		surface[i].x = q[ptr++];
		surface[i].y = q[ptr++];
		surface[i].z = q[ptr++];
	}
}

void rbspsurf(const sm::vec3* ctl_pts, int order_u, int order_v,
	          int npts, int mpts, int p1, int p2, std::vector<sm::vec3>& surface)
{
	float b[61];
	float q[301];

	memset(b, 0, sizeof(b));
	memset(q, 0, sizeof(q));

	b[1] = -15.;
	b[2] = 0.;
	b[3] = 15.;
	b[4] = 1;
	b[5] = -15.;
	b[6] = 5.;
	b[7] = 5.;
	b[8] = 1;
	b[9] = -15.;
	b[10] = 5.;
	b[11] = -5.;
	b[12] = 1;
	b[13] = -15.;
	b[14] = 0.;
	b[15] = -15.;
	b[16] = 1;

	b[17] = -5.;
	b[18] = 5.;
	b[19] = 15.;
	b[20] = 1;
	b[21] = -5.;
	b[22] = 10.;
	b[23] = 5.;
	b[24] = 1;
	b[25] = -5.;
	b[26] = 10.;
	b[27] = -5.;
	b[28] = 1;
	b[29] = -5.;
	b[30] = 5.;
	b[31] = -15.;
	b[32] = 1;

	b[33] = 5.;
	b[34] = 5.;
	b[35] = 15.;
	b[36] = 1;
	b[37] = 5.;
	b[38] = 10.;
	b[39] = 5.;
	b[40] = 1;
	b[41] = 5.;
	b[42] = 10.;
	b[43] = -5.;
	b[44] = 1;
	b[45] = 5.;
	b[46] = 0.;
	b[47] = -15.;
	b[48] = 1;

	b[49] = 15.;
	b[50] = 0.;
	b[51] = 15.;
	b[52] = 1;
	b[53] = 15.;
	b[54] = 5.;
	b[55] = 5.;
	b[56] = 1;
	b[57] = 15.;
	b[58] = 5.;
	b[59] = -5.;
	b[60] = 1;
	b[61] = 15.;
	b[62] = 0.;
	b[63] = -15.;
	b[64] = 1;

	for (int i = 0; i < 64; ++i) {
		b[i] = b[i + 1] * 0.01f;;
	}

	int n = npts - 1;
	int m = mpts - 1;
	//surface.resize(n * p1 * m * p2);
	//aitn::bezsurf(const_cast<float*>(ctl_pts[0].xyz), n, m, 
	//	p1, p2, surface[0].xyz);

	aitn::rbspsurf(b, order_u, order_v, n, m, p1, p2, q);

	surface.resize(p1 * p2);

	int ptr = 0;
	for (int i = 0, n = p1 * p2; i < n; ++i) 
	{
		surface[i].x = q[ptr++];
		surface[i].y = q[ptr++];
		surface[i].z = q[ptr++];
	}
}

}