`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date:    15:45:55 01/09/2018 
// Design Name: 
// Module Name:    add_16 
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
module add_16(
    input [15:0] a,
    input [15:0] b,
    output [15:0] s
    );
wire c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16;

FA f1 (a[0],b[0],1'b0,s[0],c1);
FA f2 (a[1],b[1],c1,s[1],c2);
FA f3 (a[2],b[2],c2,s[2],c3);
FA f4 (a[3],b[3],c3,s[3],c4);
FA f5 (a[4],b[4],c4,s[4],c5);
FA f6 (a[5],b[5],c5,s[5],c6);
FA f7 (a[6],b[6],c6,s[6],c7);
FA f8 (a[7],b[7],c7,s[7],c8);
FA f9 (a[8],b[8],c8,s[8],c9);
FA f10 (a[9],b[9],c9,s[9],c10);
FA f11 (a[10],b[10],c10,s[10],c11);
FA f12 (a[11],b[11],c11,s[11],c12);
FA f13 (a[12],b[12],c12,s[12],c13);
FA f14 (a[13],b[13],c13,s[13],c14);
FA f15 (a[14],b[14],c14,s[14],c15);
FA f16 (a[15],b[15],c15,s[15],c16);

endmodule
