/*
 * This file is part of the MicroPython project, http://micropython.org/
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2019 Dan Halbert for Adafruit Industries
 * Copyright (c) 2018 Artur Pacholec
 * Copyright (c) 2016 Glenn Ruben Bakke
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

#include "shared-bindings/bleio/__init__.h"
#include "shared-bindings/bleio/Address.h"
#include "shared-bindings/bleio/Attribute.h"
#include "shared-bindings/bleio/Central.h"
#include "shared-bindings/bleio/Characteristic.h"
#include "shared-bindings/bleio/CharacteristicBuffer.h"
#include "shared-bindings/bleio/Descriptor.h"
#include "shared-bindings/bleio/Peripheral.h"
#include "shared-bindings/bleio/ScanEntry.h"
#include "shared-bindings/bleio/Scanner.h"
#include "shared-bindings/bleio/Service.h"
#include "shared-bindings/bleio/UUID.h"

//| :mod:`bleio` --- Bluetooth Low Energy (BLE) communication
//| ================================================================
//|
//| .. module:: bleio
//|   :synopsis: Bluetooth Low Energy functionality
//|   :platform: nRF
//|
//| The `bleio` module provides necessary low-level functionality for communicating
//| using Bluetooth Low Energy (BLE). We recommend you use `bleio` in conjunction
//| with the `adafruit_ble <https://circuitpython.readthedocs.io/projects/ble/en/latest/>`_
//| CircuitPython library, which builds on `bleio`, and
//| provides higher-level convenience functionality, including predefined beacons, clients,
//| servers.
//|
//| Libraries
//|
//| .. toctree::
//|     :maxdepth: 3
//|
//|     Address
//|     Adapter
//|     Attribute
//|     Central
//|     Characteristic
//|     CharacteristicBuffer
//|     Descriptor
//|     Peripheral
//|     ScanEntry
//|     Scanner
//|     Service
//|     UUID
//|
//| .. attribute:: adapter
//|
//|   BLE Adapter information, such as enabled state as well as MAC
//|   address.
//|   This object is the sole instance of `bleio.Adapter`.
//|

STATIC const mp_rom_map_elem_t bleio_module_globals_table[] = {
    { MP_ROM_QSTR(MP_QSTR___name__),             MP_ROM_QSTR(MP_QSTR_bleio) },
    { MP_ROM_QSTR(MP_QSTR_Address),              MP_ROM_PTR(&bleio_address_type) },
    { MP_ROM_QSTR(MP_QSTR_Attribute),            MP_ROM_PTR(&bleio_attribute_type) },
    { MP_ROM_QSTR(MP_QSTR_Central),              MP_ROM_PTR(&bleio_central_type) },
    { MP_ROM_QSTR(MP_QSTR_Characteristic),       MP_ROM_PTR(&bleio_characteristic_type) },
    { MP_ROM_QSTR(MP_QSTR_CharacteristicBuffer), MP_ROM_PTR(&bleio_characteristic_buffer_type) },
    { MP_ROM_QSTR(MP_QSTR_Descriptor),           MP_ROM_PTR(&bleio_descriptor_type) },
    { MP_ROM_QSTR(MP_QSTR_Peripheral),           MP_ROM_PTR(&bleio_peripheral_type) },
    { MP_ROM_QSTR(MP_QSTR_ScanEntry),            MP_ROM_PTR(&bleio_scanentry_type) },
    { MP_ROM_QSTR(MP_QSTR_Scanner),              MP_ROM_PTR(&bleio_scanner_type) },
    { MP_ROM_QSTR(MP_QSTR_Service),              MP_ROM_PTR(&bleio_service_type) },
    { MP_ROM_QSTR(MP_QSTR_UUID),                 MP_ROM_PTR(&bleio_uuid_type) },

    // Properties
    { MP_ROM_QSTR(MP_QSTR_adapter),              MP_ROM_PTR(&common_hal_bleio_adapter_obj) },
};

STATIC MP_DEFINE_CONST_DICT(bleio_module_globals, bleio_module_globals_table);

const mp_obj_module_t bleio_module = {
    .base = { &mp_type_module },
    .globals = (mp_obj_dict_t*)&bleio_module_globals,
};
