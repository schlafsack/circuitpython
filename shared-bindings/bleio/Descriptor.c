/*
 * This file is part of the MicroPython project, http://micropython.org/
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2019 Dan Halbert for Adafruit Industries
 * Copyright (c) 2018 Artur Pacholec
 * Copyright (c) 2017 Glenn Ruben Bakke
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "py/objproperty.h"
#include "py/runtime.h"
#include "shared-bindings/bleio/Attribute.h"
#include "shared-bindings/bleio/Descriptor.h"
#include "shared-bindings/bleio/UUID.h"

//| .. currentmodule:: bleio
//|
//| :class:`Descriptor` -- BLE descriptor
//| =========================================================
//|
//| Stores information about a BLE descriptor.
//| Descriptors are encapsulated by BLE characteristics and provide contextual
//| information about the characteristic.
//|

//| .. class:: Descriptor(uuid, *, read_perm=`Attribute.OPEN`, write_perm=`Attribute.OPEN`)
//|
//|   Create a new descriptor object with the UUID uuid
//|
//|   :param bleio.UUID uuid: The uuid of the descriptor
//|   :param int read_perm: Specifies whether the descriptor can be read by a client, and if so, which
//|      security mode is required. Must be one of the integer values `Attribute.NO_ACCESS`, `Attribute.OPEN`,
//|      `Attribute.ENCRYPT_NO_MITM`, `Attribute.ENCRYPT_WITH_MITM`, `Attribute.LESC_ENCRYPT_WITH_MITM`,
//|      `Attribute.SIGNED_NO_MITM`, or `Attribute.SIGNED_WITH_MITM`.
//|   :param int write_perm: Specifies whether the descriptor can be written by a client, and if so, which
//|      security mode is required. Values allowed are the same as ``read_perm``.
//|   :param int max_length: Maximum length in bytes of the characteristic value. The maximum allowed is
//|      is 512, or possibly 510 if ``fixed_length`` is False. The default, 20, is the maximum
//|      number of data bytes that fit in a single BLE 4.x ATT packet.
//|   :param bool fixed_length: True if the characteristic value is of fixed length.
//|
STATIC mp_obj_t bleio_descriptor_make_new(const mp_obj_type_t *type, size_t n_args, const mp_obj_t *pos_args, mp_map_t *kw_args) {
    enum { ARG_uuid, ARG_read_perm, ARG_write_perm, ARG_max_length, ARG_fixed_length };
    static const mp_arg_t allowed_args[] = {
        { MP_QSTR_uuid, MP_ARG_REQUIRED | MP_ARG_OBJ },
        { MP_QSTR_read_perm, MP_ARG_KW_ONLY| MP_ARG_INT, {.u_int = SECURITY_MODE_OPEN } },
        { MP_QSTR_write_perm, MP_ARG_KW_ONLY| MP_ARG_INT, {.u_int = SECURITY_MODE_OPEN } },
        { MP_QSTR_max_length, MP_ARG_KW_ONLY| MP_ARG_INT, {.u_int = 20} },
        { MP_QSTR_fixed_length, MP_ARG_KW_ONLY| MP_ARG_BOOL, {.u_bool = false} },
    };

    mp_arg_val_t args[MP_ARRAY_SIZE(allowed_args)];
    mp_arg_parse_all(n_args, pos_args, kw_args, MP_ARRAY_SIZE(allowed_args), allowed_args, args);

    const mp_obj_t uuid_arg = args[ARG_uuid].u_obj;

    if (!MP_OBJ_IS_TYPE(uuid_arg, &bleio_uuid_type)) {
        mp_raise_ValueError(translate("Expected a UUID"));
    }

    const bleio_attribute_security_mode_t read_perm = args[ARG_read_perm].u_int;
    common_hal_bleio_attribute_security_mode_check_valid(read_perm);

    const bleio_attribute_security_mode_t write_perm = args[ARG_write_perm].u_int;
    common_hal_bleio_attribute_security_mode_check_valid(write_perm);

    bleio_descriptor_obj_t *self = m_new_obj(bleio_descriptor_obj_t);
    self->base.type = type;
    bleio_uuid_obj_t *uuid = MP_OBJ_TO_PTR(uuid_arg);

    // Range checking on max_length arg is done by the common_hal layer, because
    // it may vary depending on underlying BLE implementation.
    common_hal_bleio_descriptor_construct(self, uuid, read_perm, write_perm,
                                          args[ARG_max_length].u_int, args[ARG_fixed_length].u_bool);

    return MP_OBJ_FROM_PTR(self);
}

//|   .. attribute:: uuid
//|
//|     The descriptor uuid. (read-only)
//|
STATIC mp_obj_t bleio_descriptor_get_uuid(mp_obj_t self_in) {
    bleio_descriptor_obj_t *self = MP_OBJ_TO_PTR(self_in);

    bleio_uuid_obj_t *uuid = common_hal_bleio_descriptor_get_uuid(self);
    return uuid ? MP_OBJ_FROM_PTR(uuid) : mp_const_none;
}
MP_DEFINE_CONST_FUN_OBJ_1(bleio_descriptor_get_uuid_obj, bleio_descriptor_get_uuid);

const mp_obj_property_t bleio_descriptor_uuid_obj = {
    .base.type = &mp_type_property,
    .proxy = {(mp_obj_t)&bleio_descriptor_get_uuid_obj,
              (mp_obj_t)&mp_const_none_obj,
              (mp_obj_t)&mp_const_none_obj},
};

//|   .. attribute:: characteristic (read-only)
//|
//|     The Characteristic this Descriptor is a part of. None if not yet assigned to a Characteristic.
//|
STATIC mp_obj_t bleio_descriptor_get_characteristic(mp_obj_t self_in) {
    bleio_descriptor_obj_t *self = MP_OBJ_TO_PTR(self_in);

    return common_hal_bleio_descriptor_get_characteristic(self);
}
STATIC MP_DEFINE_CONST_FUN_OBJ_1(bleio_descriptor_get_characteristic_obj, bleio_descriptor_get_characteristic);

const mp_obj_property_t bleio_descriptor_characteristic_obj = {
    .base.type = &mp_type_property,
    .proxy = { (mp_obj_t)&bleio_descriptor_get_characteristic_obj,
               (mp_obj_t)&mp_const_none_obj,
               (mp_obj_t)&mp_const_none_obj },
};

//|   .. attribute:: value
//|
//|     The value of this descriptor.
//|
STATIC mp_obj_t bleio_descriptor_get_value(mp_obj_t self_in) {
    bleio_descriptor_obj_t *self = MP_OBJ_TO_PTR(self_in);

    return common_hal_bleio_descriptor_get_value(self);
}
STATIC MP_DEFINE_CONST_FUN_OBJ_1(bleio_descriptor_get_value_obj, bleio_descriptor_get_value);

STATIC mp_obj_t bleio_descriptor_set_value(mp_obj_t self_in, mp_obj_t value_in) {
    bleio_descriptor_obj_t *self = MP_OBJ_TO_PTR(self_in);

    mp_buffer_info_t bufinfo;
    mp_get_buffer_raise(value_in, &bufinfo, MP_BUFFER_READ);

    common_hal_bleio_descriptor_set_value(self, &bufinfo);

    return mp_const_none;
}
STATIC MP_DEFINE_CONST_FUN_OBJ_2(bleio_descriptor_set_value_obj, bleio_descriptor_set_value);

const mp_obj_property_t bleio_descriptor_value_obj = {
    .base.type = &mp_type_property,
    .proxy = { (mp_obj_t)&bleio_descriptor_get_value_obj,
               (mp_obj_t)&bleio_descriptor_set_value_obj,
               (mp_obj_t)&mp_const_none_obj },
};

STATIC const mp_rom_map_elem_t bleio_descriptor_locals_dict_table[] = {
    // Properties
    { MP_ROM_QSTR(MP_QSTR_uuid), MP_ROM_PTR(&bleio_descriptor_uuid_obj) },
    { MP_ROM_QSTR(MP_QSTR_characteristic), MP_ROM_PTR(&bleio_descriptor_characteristic_obj) },
    { MP_ROM_QSTR(MP_QSTR_value), MP_ROM_PTR(&bleio_descriptor_value_obj) },
};

STATIC MP_DEFINE_CONST_DICT(bleio_descriptor_locals_dict, bleio_descriptor_locals_dict_table);

STATIC void bleio_descriptor_print(const mp_print_t *print, mp_obj_t self_in, mp_print_kind_t kind) {
    bleio_descriptor_obj_t *self = MP_OBJ_TO_PTR(self_in);
    if (self->uuid) {
        mp_printf(print, "Descriptor(");
        bleio_uuid_print(print, MP_OBJ_FROM_PTR(self->uuid), kind);
        mp_printf(print, ")");
    } else {
        mp_printf(print, "<Descriptor with Unregistered UUID>");
    }
}

const mp_obj_type_t bleio_descriptor_type = {
    { &mp_type_type },
    .name = MP_QSTR_Descriptor,
    .make_new = bleio_descriptor_make_new,
    .print = bleio_descriptor_print,
    .locals_dict = (mp_obj_dict_t*)&bleio_descriptor_locals_dict
};
