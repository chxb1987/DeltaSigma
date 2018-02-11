`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    10:02:08 12/18/2017 
// Design Name: 
// Module Name:    diff 
// Project Name: 
// Target Devices: 
// Tool versions: 
// Description: 
//
// Dependencies: 
//
// Revision: 
// Revision 0.01 - File Created
// Additional Comments: 
//
//////////////////////////////////////////////////////////////////////////////////
module diff(
    input din,
	 input reset,
    input clk,
    output reg [11:0] dout
    );
reg [11:0] ip;
reg [11:0] diff1;
reg [11:0] ip_1;
reg [11:0] diff1_1;
always @(din)
begin
	if (din == 0)
		begin
		ip <= 0;
		end
	else
		begin
		ip <= 1;
		end
end

always @(posedge clk or posedge reset)
begin
	if (reset)
		begin
		diff1 <= 0;
		dout <= 0;
		ip_1 <= 0;
		diff1_1 <= 0;
		end
	else
		begin
		diff1 <= ip - ip_1;
		dout <= diff1 - diff1_1;
		ip_1 <= ip;
		diff1_1 <= diff1;	
		end
end

endmodule
