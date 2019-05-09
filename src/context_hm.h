#pragma once
// *******************************************************************************************
// This file is a part of GTShark software distributed under GNU GPL 3 licence.
// The homepage of the GTShark project is https://github.com/refresh-bio/GTShark
//
// Author : Sebastian Deorowicz and Agnieszka Danek
// Version: 1.1
// Date   : 2019-05-09
// *******************************************************************************************

#include <mmintrin.h>
#include <cstdint>
#include <xmmintrin.h>
#include <iostream> 
#include <cstddef>

#include "defs.h"
#include "rc.h"
#include "io.h"

template<typename MODEL> class CContextHM {
public:
	typedef struct {
		context_t ctx;
		MODEL *rcm;
	} item_t;

	typedef context_t key_type;
	typedef MODEL* value_type;

private:
	double max_fill_factor;

	size_t size;
	size_t filled;
	item_t *data;
	size_t allocated;
	size_t size_when_restruct;
	size_t allocated_mask;

	size_t ht_memory;
	size_t ht_total;
	size_t ht_match;

	void restruct(void)
	{
		item_t *old_data = data;
		size_t old_allocated = allocated;

		allocated *= 2;
		size = 0;
		filled = 0;

		allocated_mask = allocated - 1ull;
		size_when_restruct = (size_t)(allocated * max_fill_factor);

		data = new item_t[allocated];
		for (size_t i = 0; i < allocated; ++i)
			data[i].rcm = nullptr;

		ht_memory += allocated * sizeof(item_t);

		for (size_t i = 0; i < old_allocated; ++i)
			if (old_data[i].rcm != nullptr)
				insert(old_data[i].ctx, old_data[i].rcm);

		delete[] old_data;
		ht_memory -= old_allocated * sizeof(item_t);
	}

	size_t hash(context_t ctx)
	{
		return (0x9e3779b97f4a7c13ull * ctx) & allocated_mask;
	}

public:
	CContextHM()
	{
		ht_memory = 0;
		ht_total = 0;
		ht_match = 0;

		allocated = 16;
		allocated_mask = allocated - 1;

		size = 0;
		filled = 0;
		data = new item_t[allocated];
		for (size_t i = 0; i < allocated; ++i)
			data[i].rcm = nullptr;

		max_fill_factor = 0.4;

		ht_memory += allocated * sizeof(item_t);

		size_when_restruct = (size_t)(allocated * max_fill_factor);
	}

	~CContextHM()
	{
		if (data == nullptr)
			return;

		for (size_t i = 0; i < allocated; ++i)
			if (data[i].rcm)
				delete data[i].rcm;
		delete[] data;
	}

	size_t get_bytes() const {
		return ht_memory;
	}

	// Mozna to przyspieszyc tak, zebyinsert wykorzystywal wiedze o tym gdzie skonczyl szukac find
	bool insert(context_t ctx, MODEL *rcm)
	{
		if (size >= size_when_restruct)
			restruct();

		size_t h = hash(ctx);

		if (data[h].rcm != nullptr)
		{
			do
			{
				h = (h + 1) & allocated_mask;
			} while (data[h].rcm != nullptr);
		}

		if (data[h].rcm == nullptr)
			++size;

		++filled;

		data[h].ctx = ctx;
		data[h].rcm = rcm;

		return true;
	}

	MODEL* find(context_t ctx) 
	{
		size_t h = hash(ctx);

		if (data[h].rcm == nullptr)
			return nullptr;

		if (data[h].ctx == ctx)
			return data[h].rcm;

		h = (h + 1) & allocated_mask;

		while (data[h].rcm != nullptr)
		{
			if (data[h].ctx == ctx)
				return data[h].rcm;
			h = (h + 1) & allocated_mask;
		}

		return nullptr;
	}

	void prefetch(context_t ctx)
	{
		size_t h = hash(ctx);

#ifdef _WIN32
		_mm_prefetch((const char*)(data + h), _MM_HINT_T0);
#else
		__builtin_prefetch(data + h);
#endif
	}

	size_t get_size(void) const
	{
		return filled;
	}
}; 

// EOF
