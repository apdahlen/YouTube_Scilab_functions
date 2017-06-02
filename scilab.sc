// These Scilab helper functions are used for phasor manipulation to solve electrical
// engineering problems. They are especially useful for electromechanical problems
// involving transformers and motors.
//
// Notes:
//
// 1) Written for Scilab. Why Scilab? Why not? It is a common scientific language that
//    will run on Linux, Windows and Mac. It's an easy to use scripting language with
//    command line like interface. Also, the price is right...
//
// 2) We could argue that phasors values are assumed to be RMS values. In fact, that
//    is what I have done in this code. All internal representations are in RMS form.
//    Use these functions with caution as many textbook authors use peak voltages
//    especially in the introductory sections before they have introduced power. This
//    statement is is true for introductory Electrical Engineering with texts such as
//    Alexander and Irwin.
//
// 3) All operations could have been carried out using a hand held calculator such
//    as my trusty HP 35s.
//
// 4) Scilab trig functions operate using radians. However these function are written
//    using degrees. Helping functions deeg2rad() and rad2deg() are included.
//


// Electrical engineering constants

	e = %e
	j = %i



// Radian and degree conversions

	function rad = deg2rad(deg)
		rad = (%pi*deg) / 180
	endfunction

	function deg = rad2deg(rad)
		deg = (rad / %pi) * 180
	endfunction




// RMS and peak conversions

	function out = RMS2peak(in)
		out = in * 2^0.5
	endfunction

	function out = peak2RMS(in)
		out = in / 2^0.5
	endfunction



// Functions used to enter sources

	function x = source_RMS_rect(in_real, in_imag)
		x = complex(in_real, in_imag)
	endfunction

	function x = source_peak_rect(in_real, in_imag)
		x = complex(in_real, in_imag) / (2^0.5)
	endfunction

	function x = source_RMS_polar(mag, deg)
		rad = deg2rad(deg)
		unit_vector = complex(cos(rad), sin(rad))
		x = mag * unit_vector
	endfunction

	function x = source_peak_polar(mag, deg)
		rad = deg2rad(deg)
		unit_vector = complex(cos(rad), sin(rad))
		x = (mag * unit_vector) / (2^0.5)
	endfunction




// Power factor to angle conversion

	function angle_deg = pf2deg(x)
		if((x<0) | (x>1))
			error(" Error: Invalid power factor")
		end
		angle_deg = rad2deg(acos(x))
		printf("\n\t*** Caution *** set result to neg if capacitive\n")
	endfunction

	function power_factor = deg2pf(x)
		power_factor = cos(deg2rad(x))
	endfunction




// Convert reactive component to impence and vice versa

	function impence = calc_cap_Z(value, frequency)
		omega = 2 * %pi * frequency
		X_C = 1 / (omega * value)
		impence = complex(0, -X_C)
	endfunction

	function impence = calc_ind_Z(value, frequency)
		omega = 2 * %pi * frequency
		X_L = omega * value
		impence = complex( 0,  X_L)
	endfunction

	function F = calc_cap_value(Z, frequency)
		omega = 2 * %pi * frequency
		X_C = -imag(Z)
		F = 1 /(omega * X_C)
	endfunction

	function H = calc_ind_value(Z, frequency)
		omega = 2 * %pi * frequency
		X_L = imag(Z)
		H = X_L / omega
	endfunction


// These function are used to display results in an easy to understand format

	function display_polar(x)
		mag = abs(x)
		rad = atan(imag(x), real(x))
		deg = rad2deg(rad)
		printf ("    %0.2f ∠ %0.2f", mag, deg)
	endfunction

	function display_current(x, frequency)
		mag = abs(x)
		rad = atan(imag(x), real(x))
		deg = rad2deg(rad)
		printf ("    %0.2f + %0.2fj A_RMS\n",real(x), imag(x))
		printf ("    %0.2f ∠ %0.2f° A_RMS\n", mag, deg)
		printf ("    %0.2f ∠ %0.2f° A_peak\n", RMS2peak(mag), deg)
		printf ("    %0.2fe^(j%0.2f°) A_peak \n", RMS2peak(mag), deg)
		printf ("    %0.2f[COS(%0.2f°) + jSIN(%0.2f°)] A_peak\n", RMS2peak(mag), deg, deg)
		printf ("    V(t) = %0.2f COS(2π%st + %0.2f°) A_peak ", RMS2peak(mag), string(frequency), deg)
	endfunction

	function display_voltage(x, frequency)
		mag = abs(x)
		rad = atan(imag(x), real(x))
		deg = rad2deg(rad)
		printf ("    %0.2f + %0.2fj V_RMS\n",real(x), imag(x))
		printf ("    %0.2f ∠ %0.2f° V_RMS\n", mag, deg)
		printf ("    %0.2f ∠ %0.2f° V_peak\n", RMS2peak(mag), deg)
		printf ("    %0.2fe^(j%0.2f°) V_peak \n", RMS2peak(mag), deg)
		printf ("    %0.2f[COS(%0.2f°) + jSIN(%0.2f°)] V_peak\n", RMS2peak(mag), deg, deg)
		printf ("    V(t) = %0.2f COS(2π%st + %0.2f°) V_peak", RMS2peak(mag), string(frequency), deg)
	endfunction
	
	function display_complex_power(x)
		printf("\tReal power (P) = %0.2f W\n", real(x))
		printf("\tReactive power (Q) = %0.2f VAR \n", imag(x))
		
		printf ("\tComplex power (S) = %0.2f + %0.2fj VA\n",real(x), imag(x))
		mag = abs(x)
		rad = atan(imag(x), real(x))
		deg = rad2deg(rad)
		printf ("\tComplex power (S) = %0.2f ∠ %0.2f VA\n", mag, deg)
		if (deg == 0)
			printf ("\Resistive with a power factor of 1\n")
		elseif (deg > 0)
			printf ("\tPower factor = %0.2f (Inductive)\n", deg2pf(deg))
		else
			printf ("\tPower factor = %0.2f (Capacitive) \n", deg2pf(deg))
		end
	endfunction






// Vector plotting function

	function prep_plot

		a = get("current_axes")				//get the handle of the newly created axes
		a.data_bounds=[-200,-200;200,200]; //set the boundary values for the x, 
		a.axes_visible="on"; 				// makes the axes visible
		a.y_label.text = "imag"
		a.X_label.text = "real"

		//champ(real(start), imag(start), real(vector), imag(vector))
		
	//	xsegs([real(start), real(vector)] , [imag(start), imag(vector)])

		c=a.children
		c.thickness=2;
		xgrid

	endfunction


	function plot_vector(start, vector)
		a = get("current_axes")				//get the handle of the newly created axes
		xsegs([real(start), real(vector)] , [imag(start), imag(vector)])
		a.data_bounds=[-200,-200;200,200]; //set the boundary values for the x,
		c=a.children
		c.thickness=2;
	endfunction



























clf();
t=[real(start)	1     1
   -1	4     5    ];

plot(t) // plots each t column versus row size











angle_2_pf(x)



	function plot_vector(start, vector)

	a = get("current_axes")				//get the handle of the newly created axes
	a.data_bounds=[-200,-200,-200,200]; //set the boundary values for the x, 
	a.axes_visible="on"; 				// makes the axes visible
	a.y_label.text = "imag"
	a.X_label.text = "real"

	//champ(real(start), imag(start), real(vector), imag(vector))
	
	//	xsegs([real(start), real(vector)] , [imag(start), imag(vector)])

	c=a.children
	c.thickness=2;
	xgrid

endfunction





//rads = atan(imag(x), real(x))
