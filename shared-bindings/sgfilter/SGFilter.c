#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include "lib/utils/context_manager_helpers.h"
#include "py/objproperty.h"
#include "py/runtime.h"
#include "py/runtime0.h"
#include "shared-bindings/sgfilter/SGFilter.h"
#include "shared-bindings/util.h"
#include "lib/utils/buffer_helper.h"

//| .. currentmodule:: sgfilter
//|
//| :class:`SGFilter` -- A Savitzky-Golay Digital Filter
//| ====================================================================================
//|
//| Provides a smoothing filter for digital data
//|
//| .. class:: SGFilter()
//|
//|   Create an object.

STATIC mp_obj_t sgfilter_sgfilter_make_new(const mp_obj_type_t *type, size_t n_args, const mp_obj_t *pos_args, mp_map_t *kw_args) {
    //mp_int_t nl, mp_int_t nr, mp_int_t ld, mp_int_t m
    enum { ARG_nl, ARG_nr, ARG_ld, ARG_m };
    static const mp_arg_t allowed_args[] = {
        { MP_QSTR_nl, MP_ARG_INT | MP_ARG_KW_ONLY, {.u_int = 15} },
        { MP_QSTR_nr, MP_ARG_INT | MP_ARG_KW_ONLY, {.u_int = 15} },
        { MP_QSTR_ld, MP_ARG_INT | MP_ARG_KW_ONLY, {.u_int = 0} },
        { MP_QSTR_m, MP_ARG_INT | MP_ARG_KW_ONLY, {.u_int = 4} },
    };
    mp_arg_val_t args[MP_ARRAY_SIZE(allowed_args)];
    mp_arg_parse_all(n_args, pos_args, kw_args, MP_ARRAY_SIZE(allowed_args), allowed_args, args);

    mp_int_t nl = args[ARG_nl].u_int;
    mp_int_t nr = args[ARG_nr].u_int;
    mp_int_t ld = args[ARG_ld].u_int;
    mp_int_t m = args[ARG_m].u_int;

    sgfilter_sgfilter_obj_t *self = m_new_obj(sgfilter_sgfilter_obj_t);
    self->base.type = &sgfilter_sgfilter_type;
    shared_module_sgfilter_sgfilter_construct(self, nl, nr, ld, m);
    return MP_OBJ_FROM_PTR(self);
}

//|   .. method:: deinit()
//|
//|      Deinitializes the filter and releases any hardware resources for reuse.
//|
STATIC mp_obj_t sgfilter_sgfilter_deinit(mp_obj_t self_in) {
  shared_module_sgfilter_sgfilter_deinit(self_in);
  return mp_const_none;
}
STATIC MP_DEFINE_CONST_FUN_OBJ_1(sgfilter_sgfilter_deinit_obj, sgfilter_sgfilter_deinit);

STATIC void check_for_deinit(sgfilter_sgfilter_obj_t *self) {
    if (shared_module_sgfilter_sgfilter_deinited(self)) {
        raise_deinited_error();
    }
}

//|   .. method:: __enter__()
//|
//|      No-op used by Context Managers.
//|
//  Provided by context manager helper.

//|   .. method:: __exit__()
//|
//|      Automatically deinitializes the hardware when exiting a context. See
//|      :ref:`lifetime-and-contextmanagers` for more info.
//|
STATIC mp_obj_t sgfilter_sgfilter_obj___exit__(size_t n_args, const mp_obj_t *args) {
  shared_module_sgfilter_sgfilter_deinit(args[0]);
  return mp_const_none;
}
STATIC MP_DEFINE_CONST_FUN_OBJ_VAR_BETWEEN(sgfilter_sgfilter___exit___obj, 4, 4, sgfilter_sgfilter_obj___exit__);

//|   .. attribute:: filter
//|
//|     Smooths the data provided
//|
STATIC mp_obj_t sgfilter_sgfilter_obj_filter(size_t n_args, const mp_obj_t *pos_args, mp_map_t *kw_args) {

  mp_check_self(MP_OBJ_IS_TYPE(pos_args[0], &sgfilter_sgfilter_type));
  sgfilter_sgfilter_obj_t *self = MP_OBJ_TO_PTR(pos_args[0]);
  check_for_deinit(self);

  enum { ARG_data, ARG_start, ARG_end};
  static const mp_arg_t allowed_args[] = {
    { MP_QSTR_data, MP_ARG_REQUIRED | MP_ARG_OBJ },
  };

  mp_arg_val_t args[MP_ARRAY_SIZE(allowed_args)];
  mp_arg_parse_all(n_args - 1, pos_args + 1, kw_args, MP_ARRAY_SIZE(allowed_args), allowed_args, args);

  if (!MP_OBJ_IS_TYPE(args[ARG_data].u_obj, &mp_type_list)
    && !MP_OBJ_IS_TYPE(args[ARG_data].u_obj, &mp_type_tuple)
    && args[ARG_data].u_obj != mp_const_none) {
    mp_raise_ValueError(translate("data must be a list, tuple, deque or None"));
  }

  // Fetch the input data as an array of objects
  size_t data_len;
  mp_obj_t *data;
  mp_obj_get_array(args[ARG_data].u_obj, &data_len, &data);

  if(data_len < (size_t)(self->nr + self->nr + 1)) {
    mp_raise_ValueError(translate("the number of data points must be at least the size of the smoothing window."));
  }

  // Batshit SG lib uses 1 based indexes!
  mp_float_t * input = m_new(mp_float_t, data_len + 2);
  mp_float_t * output = m_new(mp_float_t, data_len + 2);

  for (size_t k = 1; k <= data_len; k++) input[k] = mp_obj_get_float(data[k-1]);

  shared_module_sgfilter_sgfilter_filter(self, input, output, data_len);

  mp_obj_t res = mp_obj_new_list(0, NULL);
  for (size_t k = 1; k <= data_len; k++) mp_obj_list_append(res, mp_obj_new_float(output[k]));

  m_del(m_float_t, input, data_len + 2);
  m_del(m_float_t, output, data_len + 2);

  return res;
}
MP_DEFINE_CONST_FUN_OBJ_KW(sgfilter_sgfilter_filter_obj, 1, sgfilter_sgfilter_obj_filter);

STATIC const mp_rom_map_elem_t sgfilter_sgfilter_locals_dict_table[] = {
    // Methods
    { MP_ROM_QSTR(MP_QSTR_deinit), MP_ROM_PTR(&sgfilter_sgfilter_deinit_obj) },
    { MP_ROM_QSTR(MP_QSTR___enter__), MP_ROM_PTR(&default___enter___obj) },
    { MP_ROM_QSTR(MP_QSTR___exit__), MP_ROM_PTR(&sgfilter_sgfilter___exit___obj) },
    { MP_ROM_QSTR(MP_QSTR_filter), MP_ROM_PTR(&sgfilter_sgfilter_filter_obj) },
};
STATIC MP_DEFINE_CONST_DICT(sgfilter_sgfilter_locals_dict, sgfilter_sgfilter_locals_dict_table);

const mp_obj_type_t sgfilter_sgfilter_type = {
    { &mp_type_type },
    .name = MP_QSTR_SGFilter,
    .make_new = sgfilter_sgfilter_make_new,
    .locals_dict = (mp_obj_dict_t*)&sgfilter_sgfilter_locals_dict,
};