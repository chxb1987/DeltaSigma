`timescale 1ns / 1ps

////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer:
//
// Create Date:   10:10:37 12/18/2017
// Design Name:   diff
// Module Name:   D:/LEARN/BTP/diff_test.v
// Project Name:  BTP
// Target Device:  
// Tool versions:  
// Description: 
//
// Verilog Test Fixture created by ISE for module: diff
//
// Dependencies:
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
////////////////////////////////////////////////////////////////////////////////

module diff_test;

	// Inputs
	reg din;
	reg reset;
	reg clk;

	// Outputs
	wire [11:0] dout;

	// Instantiate the Unit Under Test (UUT)
	diff uut (
		.din(din), 
		.reset(reset), 
		.clk(clk), 
		.dout(dout)
	);

	initial begin
		// Initialize Inputs
		din = 0;
		reset = 1;
		clk = 0;

		// Wait 100 ns for global reset to finish
		#20;
      reset = 0;
		din = 1;
		#60;
		din = 0;
		// Add stimulus here

	end
always begin #10; clk <= ~clk;  end    
endmodule

