#ifndef MICROPY_INCLUDED_SHARED_BINDINGS_SGFILTER_SGFILTER_H
#define MICROPY_INCLUDED_SHARED_BINDINGS_SGFILTER_SGFILTER_H

#include "shared-module/sgfilter/SGFilter.h"

extern const mp_obj_type_t sgfilter_sgfilter_type;

extern void shared_module_sgfilter_sgfilter_construct(sgfilter_sgfilter_obj_t* self, mp_int_t nl, mp_int_t nr, mp_int_t ld, mp_int_t m);
extern void shared_module_sgfilter_sgfilter_deinit(sgfilter_sgfilter_obj_t* self);
extern bool shared_module_sgfilter_sgfilter_deinited(sgfilter_sgfilter_obj_t* self);
extern void shared_module_sgfilter_sgfilter_filter(sgfilter_sgfilter_obj_t* self, mp_float_t in[], mp_float_t out[], size_t len);

#endif // MICROPY_INCLUDED_SHARED_BINDINGS_SGFILTER_SGFILTER_H