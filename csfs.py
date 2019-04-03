#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

    A small script that calculates the number of 
    configuration state functions (CSFs), e.g. 
    for CASSCF calculations

    written 2019 by Michael Böhme
    https://github.com/micb25/chemtools

"""

import math
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

class MainWindow(Gtk.Window):
	
	paramS = Gtk.Adjustment(1.5, 0.0, 100.0, 0.5, 0, 0)
	paramN = Gtk.Adjustment(7, 1, 100, 1, 0, 0)
	paramM = Gtk.Adjustment(5, 1, 100, 1, 0, 0)
	result = Gtk.Entry()
	
	def update_values(self, widget):
		self.mult = 1 + 2 * self.paramS.get_value()
		self.npo = 1 + self.paramM.get_value()
		self.k1 = self.mult / (self.paramM.get_value() + 1)
		self.f1 = int(self.paramN.get_value()/2.0 - self.paramS.get_value())
		self.f2 = int(self.npo - (self.paramN.get_value()/2.0 - self.paramS.get_value()))
		self.f3 = int(self.paramN.get_value()/2.0 + self.paramS.get_value() + 1)
		self.f4 = int(self.npo - (self.paramN.get_value()/2.0 + self.paramS.get_value() + 1))
		if (self.f1 < 0) or (self.f2 < 0) or (self.f3 < 0) or (self.f4 < 0):
			self.result.set_text("invalid input!")
		else:
			self.k2 = math.factorial(self.npo) / ( math.factorial( self.f1 ) * math.factorial( self.f2 ) )
			self.k3 = math.factorial(self.npo) / ( math.factorial( self.f3 ) * math.factorial( self.f4 ) )
			self.result.set_text(str(int(self.k1 * self.k2 * self.k3)))
		
	
	def __init__(self):
		Gtk.Window.__init__(self, title="Number of CSFs")
		self.set_position(Gtk.WindowPosition.CENTER)		
		
		self.grid = Gtk.Grid()
		self.grid.props.margin_top = 30
		self.grid.props.margin_bottom = 30
		self.grid.props.margin_left = 30
		self.grid.props.margin_right = 30
		self.grid.set_row_spacing(30)
		self.grid.set_column_spacing(30)
		self.add(self.grid)
				
		self.label1 = Gtk.Label()
		self.label1.set_markup("Spin <i>S</i> =")
		self.label1.set_justify(Gtk.Justification.RIGHT)
		self.grid.attach(self.label1, 0, 0, 1, 1)
		self.spin1 = Gtk.SpinButton(adjustment=self.paramS, climb_rate=0.5, digits=1)
		self.spin1.connect("value-changed", self.update_values )
		self.grid.attach(self.spin1, 1, 0, 1, 1)
		
		self.label2 = Gtk.Label()
		self.label2.set_markup("Number of electrons <i>n</i> =")
		self.label2.set_justify(Gtk.Justification.RIGHT)
		self.grid.attach(self.label2, 0, 1, 1, 1)
		self.spin2 = Gtk.SpinButton(adjustment=self.paramN, climb_rate=1, digits=0)
		self.spin2.connect("value-changed", self.update_values )
		self.grid.attach(self.spin2, 1, 1, 1, 1)
		
		self.label3 = Gtk.Label()
		self.label3.set_markup("Number of orbitals <i>m</i> =")
		self.label3.set_justify(Gtk.Justification.RIGHT)
		self.grid.attach(self.label3, 0, 2, 1, 1)
		self.spin3 = Gtk.SpinButton(adjustment=self.paramM, climb_rate=1, digits=0)
		self.spin3.connect("value-changed", self.update_values )
		self.grid.attach(self.spin3, 1, 2, 1, 1)
		
		self.hsep = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
		self.grid.attach(self.hsep, 1, 3, 1, 1)
		
		self.label4 = Gtk.Label()
		self.label4.set_markup("Number of CSFs =")
		self.label4.set_justify(Gtk.Justification.RIGHT)
		self.grid.attach(self.label4, 0, 4, 1, 1)
		self.result = Gtk.Entry()
		self.result.set_text("0")
		self.result.set_editable(False)
		self.grid.attach(self.result, 1, 4, 1, 1)
		
		self.update_values(self)

w = MainWindow()
w.connect("destroy", Gtk.main_quit)
w.show_all()
Gtk.main()
