#include "volume4d.h"
#include <fstream>
#include <cassert>
#include <iostream>
#include <type_traits>

template <class Ty>
Volume4d<Ty>::Volume4d(size_t nx, size_t ny, size_t nz, size_t nt) noexcept :
	data(),
	nx(nx),
	ny(ny),
	nz(nz),
	nt(nt)
{}

template <class Ty>
Volume4d<Ty>::Volume4d(const Volume4d& other) noexcept :
	data(other.data),
	nx(other.nx),
	ny(other.ny),
	nz(other.nz),
	nt(other.nt)
{}

template <class Ty>
Volume4d<Ty>::Volume4d(Volume4d&& other) noexcept :
	data(std::move(other.data)),
	nx(other.nx),
	ny(other.ny),
	nz(other.nz),
	nt(other.nt)
{
	other.nx = 0;
	other.ny = 0;
	other.nz = 0;
	other.nt = 0;
}

template <class Ty>
Volume4d<Ty>& Volume4d<Ty>::operator=(const Volume4d& other) noexcept
{
	if (this != &other) {
		data = other.data;
		nx = other.nx;
		ny = other.ny;
		nz = other.nz;
		nt = other.nt;
	}
	return *this;
}

template <class Ty>
Volume4d<Ty>& Volume4d<Ty>::operator=(Volume4d&& other) noexcept
{
	if (this != &other) {
		data = std::move(other.data);
		nx = other.nx;
		ny = other.ny;
		nz = other.nz;
		nt = other.nt;

		other.nx = 0;
		other.ny = 0;
		other.nz = 0;
		other.nt = 0;
	}
	return *this;
}

template <class Ty>
void Volume4d<Ty>::alloc()
{
	alloc(nx, ny, nz, nt);
}

template <class Ty>
void Volume4d<Ty>::alloc(size_t nelem)
{
	// Need to supply custom deleter since shared_ptr does not support
	// arrays prior to C++17
	data = std::shared_ptr<Ty>(new Ty[nelem], std::default_delete<Ty[]>());
}

template <class Ty>
void Volume4d<Ty>::alloc(size_t nx, size_t ny, size_t nz, size_t nt)
{
	alloc(nx * ny * nz * nt);
}

template <class Ty>
void Volume4d<Ty>::clear() {
	nx = 0;
	ny = 0;
	nz = 0;
	nt = 0;
	data = nullptr;
}

template <class Ty>
void Volume4d<Ty>::copy(const Volume4d<Ty>& src)
{
	if (this != &src) {
		nx = src.nx;
		ny = src.ny;
		nz = src.nz;
		nt = src.nt;
		alloc();
		std::memcpy(data.get(), src.get(), nx * ny * nz * sizeof(Ty));
	}
}

template <class Ty>
void Volume4d<Ty>::copy(Volume4d<Ty>&& src) noexcept
{
	if (this != &src) {
		// If src is an r-value we can just move from it instead of doing a deep copy
		*this = std::move(src);
	}
}

template <class Ty>
size_t Volume4d<Ty>::numElem() const noexcept
{
	return nx * ny * nz * nt;
}

template <class Ty>
Ty Volume4d<Ty>::interp3d(float x, float y, float z, size_t t) const
{
	float fx = x - floorf(x);
	float fy = y - floorf(y);
	float fz = z - floorf(z);
	size_t x0 = static_cast<size_t>(floorf(x));
	size_t y0 = static_cast<size_t>(floorf(y));
	size_t z0 = static_cast<size_t>(floorf(z));
	// Handle special case where query point is on boundary by just re-using first point
	size_t x1 = (x0 == nx - 1 && x <= nx - 1) ? x0 : x0 + 1;
	size_t y1 = (y0 == ny - 1 && y <= ny - 1) ? y0 : y0 + 1;
	size_t z1 = (z0 == nz - 1 && z <= nz - 1) ? z0 : z0 + 1;

	// Do interpolation
	Ty c000 = at(x0, y0, z0, t);
	Ty c001 = at(x0, y0, z1, t);
	Ty c010 = at(x0, y1, z0, t);
	Ty c011 = at(x0, y1, z1, t);
	Ty c100 = at(x1, y0, z0, t);
	Ty c101 = at(x1, y0, z1, t);
	Ty c110 = at(x1, y1, z0, t);
	Ty c111 = at(x1, y1, z1, t);
	Ty c00 = (1.0f - fx) * c000 + fx * c100;
	Ty c01 = (1.0f - fx) * c001 + fx * c101;
	Ty c10 = (1.0f - fx) * c010 + fx * c110;
	Ty c11 = (1.0f - fx) * c011 + fx * c111;
	Ty c0 = (1.0f - fy) * c00 + fy * c10;
	Ty c1 = (1.0f - fy) * c01 + fy * c11;
	return (1.0f - fz) * c0 + fz * c1;
}

template <class Ty>
Ty Volume4d<Ty>::interp3d(const Vec3f& p, size_t t) const
{
	return interp3d(p[0], p[1], p[2], t);
}

template <class Ty>
Ty Volume4d<Ty>::interp(float x, float y, float z, float t) const
{
	float ft = t - floorf(t);
	size_t t0 = static_cast<size_t>(floorf(t));
	// Handle special case where query point is on boundary by just re-using first point
	size_t t1 = (t0 == nt - 1 && t <= nt - 1) ? t0 : t0 + 1;

	// Do interpolation
	// TODO: This is not the most efficient way to do this - profile to see if it matters
	Ty c0 = interp3d(x, y, z, t0);
	Ty c1 = interp3d(x, y, z, t1);
	return (1.0f - ft) * c0 + ft * c1;
}

template <class Ty>
Ty Volume4d<Ty>::interp(const Vec4f& p) const
{
	return interp(p[0], p[1], p[2], p[3]);
}

template <class Ty>
size_t Volume4d<Ty>::idx(size_t x, size_t y, size_t z, size_t t) const noexcept
{
	return x + y * nx + z * nx * ny + t * nx * ny * nz;
}

template <class Ty>
Vec4i Volume4d<Ty>::pos(size_t idx) const noexcept
{
	return Vec4i(
		idx % nx,
		(idx / nx) % ny,
		(idx / (nx * ny)) % nz,
		idx / (nx * ny * nz)
	);
}

template <class Ty>
bool Volume4d<Ty>::contains(const Vec4f& p) const noexcept
{
	return 0 <= p[0] && 0 <= p[1] && 0 <= p[2] && 0 <= p[3] &&
		p[0] <= nx - 1 && p[1] <= ny - 1 && p[2] <= nz - 1 && p[3] <= nt - 1;
}

template <class Ty>
Ty& Volume4d<Ty>::at(size_t i)
{
	assert(i < nx * ny * nz * nt);
	return data.get()[i];
}

template <class Ty>
const Ty& Volume4d<Ty>::at(size_t i) const
{
	assert(i < nx * ny * nz * nt);
	return data.get()[i];
}

template <class Ty>
Ty& Volume4d<Ty>::at(size_t x, size_t y, size_t z, size_t t)
{
	return at(idx(x, y, z));
}

template <class Ty>
const Ty& Volume4d<Ty>::at(size_t x, size_t y, size_t z, size_t t) const
{
	return at(idx(x, y, z, t));
}

template <class Ty>
Ty& Volume4d<Ty>::at(const Vec4i& p)
{
	assert(p[0] >= 0 && p[1] >= 0 && p[2] >= 0 && p[3] >= 0);
	return at(p[0], p[1], p[2], p[3]);
}

template <class Ty>
const Ty& Volume4d<Ty>::at(const Vec4i& p) const
{
	assert(p[0] >= 0 && p[1] >= 0 && p[2] >= 0 && p[3] >= 0);
	return at(p[0], p[1], p[2], p[3]);
}

template <class Ty>
Ty& Volume4d<Ty>::operator[](size_t i)
{
	return at(i);
}

template <class Ty>
const Ty& Volume4d<Ty>::operator[](size_t i) const
{
	return at(i);
}

template <class Ty>
Ty& Volume4d<Ty>::operator[](const Vec4i &p)
{
	return at(p);
}

template <class Ty>
const Ty& Volume4d<Ty>::operator[](const Vec4i& p) const
{
	return at(p);
}