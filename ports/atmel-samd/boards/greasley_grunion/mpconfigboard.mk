USB_VID = 0x239A
USB_PID = 0x8021
USB_PRODUCT = "Greasley Grunion"
USB_MANUFACTURER = "Tom Greasley"

CHIP_VARIANT = SAMD51J19A
CHIP_FAMILY = samd51

QSPI_FLASH_FILESYSTEM = 1
EXTERNAL_FLASH_DEVICE_COUNT = 3
EXTERNAL_FLASH_DEVICES = "S25FL116K, S25FL216K, GD25Q16C"
LONGINT_IMPL = MPZ

MICROPY_PY_WIZNET5K = 5500

CIRCUITPY_ANALOGIO  = 0
CIRCUITPY_AUDIOBUSIO = 0
CIRCUITPY_AUDIOIO = 0
CIRCUITPY_AUDIOIO_COMPAT = 0
CIRCUITPY_AUDIOPWMIO = 0
CIRCUITPY_AUDIOCORE = 0
CIRCUITPY_AUDIOMIXER = 0
CIRCUITPY_AUDIOMP3 = 0
CIRCUITPY_BITBANGIO = 1
CIRCUITPY_BLEIO = 0
CIRCUITPY_BOARD = 1
CIRCUITPY_BUSIO = 1
CIRCUITPY_DIGITALIO = 1
CIRCUITPY_DISPLAYIO = 0
CIRCUITPY_FRAMEBUFFERIO = 0
CIRCUITPY_FREQUENCYIO = 0
CIRCUITPY_GAMEPAD = 0
CIRCUITPY_GAMEPADSHIFT = 0
CIRCUITPY_I2CSLAVE = 0
CIRCUITPY_MATH = 1
CIRCUITPY__EVE = 0
CIRCUITPY_MICROCONTROLLER = 1
CIRCUITPY_NEOPIXEL_WRITE = 1
CIRCUITPY_NETWORK = 1
CIRCUITPY_NVM = 1
CIRCUITPY_OS = 1
CIRCUITPY_PIXELBUF = 0
CIRCUITPY_RGBMATRIX = 0
CIRCUITPY_PULSEIO = 0
CIRCUITPY_PS2IO = 1
CIRCUITPY_RANDOM = 1
CIRCUITPY_ROTARYIO = 0
CIRCUITPY_RTC = 1
CIRCUITPY_SAMD = 1
CIRCUITPY_STAGE = 0
CIRCUITPY_STORAGE = 1
CIRCUITPY_STRUCT = 1
CIRCUITPY_SUPERVISOR = 1
CIRCUITPY_TIME = 1
CIRCUITPY_TOUCHIO_USE_NATIVE = 0
CIRCUITPY_TOUCHIO = 0
CIRCUITPY_UHEAP = 0
CIRCUITPY_USB_HID = 0
CIRCUITPY_USB_MIDI = 0
CIRCUITPY_PEW = 0
CIRCUITPY_USTACK = 0
CIRCUITPY_BITBANG_APA102 = 0
CIRCUITPY_REQUIRE_I2C_PULLUPS = 1
CIRCUITPY_SERIAL_BLE = 0
CIRCUITPY_BLE_FILE_SERVICE = 0
CIRCUITPY_SERIAL_UART = 0
CIRCUITPY_ENABLE_MPY_NATIVE = 0

CIRCUITPY_SGFILTER = 1

FROZEN_MPY_DIRS += $(TOP)/frozen/Adafruit_CircuitPython_BusDevice
FROZEN_MPY_DIRS += $(TOP)/frozen/Adafruit_CircuitPython_NeoPixel
FROZEN_MPY_DIRS += $(TOP)/frozen/Adafruit_CircuitPython_ESP32SPI
FROZEN_MPY_DIRS += $(TOP)/frozen/CircuitPython_NCD_PR33_15