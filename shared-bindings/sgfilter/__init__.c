#include <stdint.h>

#include "py/obj.h"
#include "py/runtime.h"

#include "shared-bindings/sgfilter/__init__.h"
#include "shared-bindings/sgfilter/SGFilter.h"

STATIC const mp_rom_map_elem_t sgfilter_module_globals_table[] = {
    { MP_ROM_QSTR(MP_QSTR___name__), MP_ROM_QSTR(MP_QSTR_sgfilter) },
    { MP_ROM_QSTR(MP_QSTR_SGFilter), MP_ROM_PTR(&sgfilter_sgfilter_type) },
};

STATIC MP_DEFINE_CONST_DICT(sgfilter_module_globals, sgfilter_module_globals_table);

const mp_obj_module_t sgfilter_module = {
    .base = { &mp_type_module },
    .globals = (mp_obj_dict_t*)&sgfilter_module_globals,
};