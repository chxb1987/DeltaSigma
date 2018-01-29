`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    15:28:26 01/09/2018 
// Design Name: 
// Module Name:    FA 
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
module FA(
    input A,
    input B,
    input C0,
    output S,
    output C1
    );
wire ABx, ABn, ACn, BCn, ot;
xor x1(ABx,A,B);
xor x2(S,ABx,C0);
and a1(ABn,A,B);
and a2(ACn,A,C0);
and a3(BCn,B,C0);
or o1(ot,ABn,ACn);
or o2(C1,ot,BCn);
endmodule
