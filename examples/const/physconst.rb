#!/usr/bin/env ruby
require("gsl")
require("curses")

#SUBWIN_H = 20
SUBWIN_W = 50
SUBWIN_Y = 1
SUBWIN_X = 10

print("MKSA/CGSM ? [m/c] : ")
sys = gets
if /^m/ =~ sys
  include GSL::CONST::MKSA
  PREFIX = "MKSA"
else
  include GSL::CONST::CGSM
  PREFIX = "CGSM"
end
include GSL::CONST::NUM

Curses.init_screen
W = Curses.stdscr
Curses.refresh

MAIN_MENU = ["", "Fundamental Constants", "Astronomy and Astrophysics",
             "Atomic and Nuclear Physics", "Measurement of Time",
             "Imperial Units", "Nautical Units", "Printers Units",
             "Volume", "Mass and Weight", "Thermal Energy and Power",
             "Pressure", "Viscosity", "Light and Illumination",
             "Radioactivity", "Force and Energy", "Prefixes"]

MENU = Array.new(MAIN_MENU.size)
MENU[1] = ["", "SPEED_OF_LIGHT", "VACUUM_PERMEABILITY", "VACUUM_PERMITTIVITY",
           "PLANCKS_CONSTANT_H", "PLANCKS_CONSTANT_HBAR", "AVOGADRO",
           "FARADAY", "BOLTZMANN", "MOLAR_GAS", "STANDARD_GAS_VOLUME",
           "STEFAN_BOLTZMANN_CONSTANT", "GAUSS", "MICRON", "HECTARE",
           "MILES_PER_HOUR", "KILOMETERS_PER_HOUR"]

MENU[2] = ["", "ASTRONOMICAL_UNIT", "GRAVITATIONAL_CONSTANT",
           "LIGHT_YEAR", "PARSEC", "GRAV_ACCEL", "SOLAR_MASS"]

MENU[3] = ["", "ELECTRON_CHARGE", "ELECTRON_VOLT", "UNIFIED_ATOMIC_MASS",
           "MASS_ELECTRON", "MASS_MUON", "MASS_PROTON", "MASS_NEUTRON",
           "FINE_STRUCTURE", "RYDBERG", "BOHR_RADIUS", "ANGSTROM",
           "BARN", "BOHR_MAGNETON", "NUCLEAR_MAGNETON",
           "ELECTRON_MAGNETIC_MOMENT", "PROTON_MAGNETIC_MOMENT",
           "THOMSON_CROSS_SECTION"]

MENU[4] = ["", "MINUTE", "HOUR", "DAY", "WEEK"]
MENU[5] = ["", "INCH", "FOOT", "YARD", "MILE", "MIL"]
MENU[6] = ["", "NAUTICAL_MILE", "FATHOM", "KNOT"]
MENU[7] = ["", "POINT", "TEXPOINT"]
MENU[8] = ["", "ACRE", "LITER", "US_GALLON", "CANADIAN_GALLON", "UK_GALLON",
           "QUART", "PINT"]
MENU[9] = ["", "POUND_MASS", "OUNCE_MASS", "TON", "METRIC_TON", "UK_TON",
           "TROY_OUNCE", "CARAT", "GRAM_FORCE", "POUND_FORCE",
           "KILOPOUND_FORCE", "POUNDAL"]
MENU[10] = ["", "CALORIE", "BTU", "THERM", "HORSEPOWER"]
MENU[11] = ["", "BAR", "STD_ATMOSPHERE", "TORR", "METER_OF_MERCURY",
            "INCH_OF_MERCURY", "INCH_OF_WATER", "PSI"]
MENU[12] = ["", "POISE", "STOKES"]
MENU[13] = ["", "STILB", "LUMEN", "LUX", "PHOT", "FOOTCANDLE",
            "LAMBERT", "FOOTLAMBERT"]
MENU[14] = ["", "CURIE", "ROENTGEN", "RAD"]
MENU[15] = ["", "NEWTON", "DYNE", "JOULE", "ERG"]
MENU[16] = ["", "YOTTA", "ZETTA", "EXA", "PETA", "TERA", "GIGA", "MEGA",
            "KILO", "MILLI", "MICRO", "NANO", "PICO", "FEMTO", "ATTO",
            "ZEPTO", "YOCTO"]

def mainmenu
  i = 1
  Curses.clear
  W.setpos(0, 1)
  W.addstr("  GSL-#{GSL::VERSION} Physical Constants:\n")
  MAIN_MENU[1..-1].each do |str|
    W.setpos(i, 1)
    W.addstr(sprintf("  [%2d] %s\n", i, str))
    i+=1
  end
  W.setpos(i, 1)
  str = "  ? ['q' to quit]: "
  W.addstr(str)
  W.setpos(i, str.size+1)
  W.refresh
end

def show_submenu(ind, w)
  w.setpos(1, 0)
  w.addstr(sprintf(" [%2d] %s\n", ind, MAIN_MENU[ind]))
  i = 1
  MENU[ind][1..-1].each do |str|
    w.setpos(i+1, 1)
    w.addstr(sprintf("    [%2d] %s\n", i, str))
    i+=1
  end
  w.setpos(i+=1, 1)
  str = "    ? ['m' to menu]: "
  w.addstr(str)
  w.setpos(i, str.size+1)
  w.box('|', '-')
  w.refresh
end

def show_result(ind, i, w)
  return if i == 0
  return if i > MENU[ind].size-1
  w.setpos(MENU[ind].size+2, 0)
  w.addstr("  #{MENU[ind][i]} = ")
  begin
    val = eval("#{MENU[ind][i]}")
    w.addstr("#{val}\n")
  rescue NameError
    w.addstr("not defined.\n")
  rescue
    return
  ensure
    w.refresh
  end
end

def submenu(ind)
  return unless MENU[ind]
  w = W.subwin(MENU[ind].size+4, SUBWIN_W, SUBWIN_Y, SUBWIN_X)
  while true
    show_submenu(ind, w)
    i = w.getstr
    if /^m/ =~ i
      w.close
      return
    end
    show_result(ind, i.to_i, w)
  end
end

#########

while true
  mainmenu
  i = W.getstr
  W.refresh
  break if /^q/ =~ i
  i = i.to_i
  next if i > MAIN_MENU.size or i < 1
  submenu(i)
end

Curses.close_screen
exit

__END__

