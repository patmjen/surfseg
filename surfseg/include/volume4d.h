#ifndef VOLUME_4D_H__
#define VOLUME_4D_H__

#include <memory>

#include <GEL/CGLA/Vec4i.h>
#include <GEL/CGLA/Vec4f.h>
#include <GEL/CGLA/Vec3i.h>

using namespace CGLA;

/** Simple structure to hold time series volumetric data of type Ty.
 * Adapted from the Volume class.
 */
template <class Ty>
struct Volume4d {
	typedef Ty DataType;

	std::shared_ptr<Ty> data;
	size_t nx;
	size_t ny;
	size_t nz;
	size_t nt;

	explicit Volume4d() :
		data(),
		nx(0),
		ny(0),
		nz(0) {}
	explicit Volume4d(size_t nx, size_t ny, size_t nz, size_t nt) noexcept;
	Volume4d(const Volume4d&) noexcept;
	Volume4d(Volume4d&& other) noexcept;

	Volume4d& operator=(const Volume4d&) noexcept;
	Volume4d& operator=(Volume4d&& other) noexcept;

	void alloc();
	void alloc(size_t nelem);
	void alloc(size_t nx, size_t ny, size_t nz, size_t nt);

	void clear();

	void copy(const Volume4d<Ty>& src);
	void copy(Volume4d<Ty>&& src) noexcept;

	size_t numElem() const noexcept;

	Ty interp3d(float x, float y, float z, size_t t) const;
	Ty interp3d(const Vec3f& p, size_t t) const;

	Ty interp(float x, float y, float z, float t) const;
	Ty interp(const Vec4f& p) const;
	
	size_t idx(size_t x, size_t y, size_t z, size_t t) const noexcept;
	Vec4i pos(size_t idx) const noexcept;
	bool contains(const Vec4f& p) const noexcept;

	Ty& at(size_t i);
	const Ty& at(size_t i) const;
	Ty& at(size_t x, size_t y, size_t z, size_t t);
	const Ty& at(size_t x, size_t y, size_t z, size_t t) const;
	Ty& at(const Vec4i& p);
	const Ty& at(const Vec4i& p) const;

	Ty& operator[](size_t i);
	const Ty& operator[](size_t i) const;
	Ty& operator[](const Vec4i& p);
	const Ty& operator[](const Vec4i& p) const;
};

#include <volume4d.inl>

#endif // VOLUME_H__
