//fir3 = FischerTechnik Industry Robots ROB3

//ECMAScript module

//npm install https://github.com/PeterTadich/singular-value-decomposition https://github.com/PeterTadich/homogeneous-transformations https://github.com/PeterTadich/matrix-computations

// To do:
//   - use 'pseudo-inverse' in IR_rob3()

import * as hlao from 'matrix-computations';
import * as svdcmp from 'singular-value-decomposition';
import * as mcht from 'homogeneous-transformations';
//import * as hlao from '../matrix-computations/hlao.mjs';
//import * as svdcmp from '../singular-value-decomposition/svdcmp.mjs';
//import * as mcht from '../homogeneous-transformations/mcht.mjs';
//import * as hlao from '../../node_modules/matrix-computations/hlao.mjs';
//import * as svdcmp from '../../node_modules/singular-value-decomposition/svdcmp.mjs';
//import * as mcht from '../../node_modules/homogeneous-transformations/mcht.mjs';

//Example:
/*
var q = [[0.0],[0.0],[0.0]]; //Home.
var RRMC = [.1,1]; //dt = 0.1, npts = 1.
var v_star = [[0.0], [-0.1], [0.0]];
var data = IR_rob3(q,RRMC,v_star); //returns [qd,T0cam]
var qd = data[0];
var T0cam = data[1];
*/
/*
pe (final):
[
    [ 0.000],
    [ 0.130],
    [-0.170]
];
*/

//   - 'q' generalised coordinates. A column vector [[v1],[d2],[d3]]
//   - 'RRMC' Resolved-Rate Motion Control - discrete time parameters. A row vector [dt,npts].
//   - 'v_star' end-effector velocity.
function IR_rob3(q,RRMC,v_star){
    //Debugger flag.
    var debug = 0;

    //Constants for the calculation of the end-effector velocity.
    //var dt = 0.1; //seconds
    //var npts = 100; //number of discrete time points.
    var dt = RRMC[0]; //seconds
    var npts = RRMC[1]; //number of discrete time points.
    
    //Define link frames using the Denavit-Hartenberg Convention.
    //var v1 = 0.0;
    //var d2 = 0.1; //100m
    //var d3 = 0.1; //100mm
    var pi =  Math.PI;
    var v1 = q[0][0]; //radians
    var d2 = 0.170; //m
    var a3 = 0.140; //m
    var d3 = q[1][0]; //m
    var v4 = -1.0*pi/2.0; //radians
    var d5 = q[2][0];
    var d6 = 0.300; //m
    var alpha_6 = -1.0*pi/2.0; //radians
    
    //Setup the links.
    var links = [
        //[ ai,alpha_i,  di,     vi]
          [0.0,      pi, 0.0,     v1], //link  1;  z0 to  z1 ( z1 = T0[0])
          [0.0,      pi,  d2, pi/2.0], //link  2;  z1 to  z2 ( z2 = T0[1])
          [ a3,     0.0,  d3,     pi], //link  3;  z2 to  z3 ( z3 = T0[2])
          [0.0,  pi/2.0, 0.0,     v4], //link  4;  z3 to  z4 ( z4 = T0[3])
          [0.0,     0.0,  d5,    0.0], //link  5;  z4 to  z5 ( z5 = T0[4])
          [0.0, alpha_6,  d6,    0.0]  //link  6;  z5 to  z6 ( z6 = T0[5])
    ];
    
    //Number of joints.
    var nbr_links = links.length;
    if(debug) console.log("Number of links: " + nbr_links);
    
    //Resolved-Rate Motion Control - discrete time 'k'.
    for(var k=0;k<npts;k=k+1){
        //Convert the Denavit-Hartenberg parameters to the homogenous transformation matrix.
        var T = [];
        T[0] = [ //z0
            [1.0,0.0,0.0,0.0],
            [0.0,1.0,0.0,0.0],
            [0.0,0.0,1.0,0.0],
            [0.0,0.0,0.0,1.0]
        ];
        for(var i=0;i<nbr_links;i=i+1){
            T.push(mcht.Aij(links[i])); //transform from from {Ai} to {Ai-1}
            if(debug){
                console.log("T[" + i + "]:");
                print_multi_array(T[i]);
            }
        }
        
        //Calculate the Jacobian matrix.
        var T0 = [];
        var z = [];
        var p = [];
        
        var z0 = [[0.0],[0.0],[1.0]]; //Selector. Equ 3.31.
        
        for(var i=0;i<=nbr_links;i=i+1){
            if(i == 0){
                T0[i] = T[i];
            } else {
                T0[i] = hlao.matrix_multiplication(T0[i-1],T[i]);
            }
            if(debug){
                console.log("T0[" + i + "]:");
                print_multi_array(T0[i]);
            }
            
            z[i] = hlao.matrix_multiplication(
                [
                    [T0[i][0][0],T0[i][0][1],T0[i][0][2]],
                    [T0[i][1][0],T0[i][1][1],T0[i][1][2]],
                    [T0[i][2][0],T0[i][2][1],T0[i][2][2]]
                ],
                z0
            );
            if(debug){
                console.log("z[" + i + "]:");
                print_multi_array(z[i]);
            }
            
            p[i] = [
                [T0[i][0][3]],
                [T0[i][1][3]],
                [T0[i][2][3]]
            ];
            if(debug){
                console.log("p[" + i + "]:");
                print_multi_array(p[i]);
            }
        }
        
        var pe = [
            [T0[nbr_links][0][3]],
            [T0[nbr_links][1][3]],
            [T0[nbr_links][2][3]]
        ];
        if(1){
            console.log("pe (initial):");
            print_multi_array(pe);
        }
        //pe = [[0],[-0.1],[0.1]];
        
        //Camera homogenesous transform with ref. to the world frame.
        var T0cam = T0[nbr_links];
        if(1){
            console.log("T0cam:");
            print_multi_array(T0cam);
        }
        
        //Calculation of the end-effector velocity.
        if((k == 0) && (typeof v_star == 'undefined')){
            //'pe' the current end-effector position (initial end-effector position) t = 0. Use forward kinematics.
            if(1){
                console.log("pe (initial end-effector position):");
                print_multi_array(pe);
            }
            //'pe_star' the desired end-effector position (final end-effector position).
            var pe_star = [[pe[0][0]+0.1],[pe[1][0]],[pe[2][0]]];
            if(1){
                console.log("pe* (final end-effector position):");
                print_multi_array(pe_star);
            }
            //End-effector velocity.
            var v_star = 
                hlao.matrix_multiplication_scalar(
                    hlao.matrix_arithmetic(pe_star,pe,'-'),
                    1.0/(dt*npts) //dt*npts = total time (dt is time step, npts is the number of time steps)
                );
            if(1){
                console.log('End-effector velocity:');
                print_multi_array(v_star);
            }
            //v_star = [[ 0.001 ],[ 0 ],[ 0 ]]
        }
        
        //Integration. page 466 (use trnorm() to ensure the transformation remains a proper homogenous transformation)
        //var Tc = trnorm(Tc * delta2tr(Vc)); //IMPORTANT: investigate this.
        
        //Transform (rotation only not translation) 'v_star' from camera coordinate frame to world coordinate frame.
        if(1){
            console.log("v_star (camera frame):");
            print_multi_array(v_star);
        }
        var Rcam = [
            [T0cam[0][0],T0cam[0][1],T0cam[0][2]],
            [T0cam[1][0],T0cam[1][1],T0cam[1][2]],
            [T0cam[2][0],T0cam[2][1],T0cam[2][2]]
        ];
        var v_star = hlao.matrix_multiplication(Rcam,v_star);
        if(1){
            console.log("v_star (world frame):");
            print_multi_array(v_star);
        }
        
        //Assemble the Jacobian matrix.
        //IMPORTANT: setup dependent on robot.
        var rj = hlao.vector_cross(z[0],hlao.matrix_arithmetic(pe,p[0],'-')); //Revolute joint.
        var JT = [ //Jacobian transpose.
            [  rj[0][0],   rj[1][0],   rj[2][0], z[0][0][0], z[0][1][0], z[0][2][0]], //zo
            [z[2][0][0], z[2][1][0], z[2][2][0],        0.0,        0.0,        0.0], //z2
            [z[4][0][0], z[4][1][0], z[4][2][0],        0.0,        0.0,        0.0]  //z4 (z5 is 'pe')
        ];
        if(debug){
            console.log("J transpose:");
            print_multi_array(JT);
        }
        
        var J = hlao.matrix_transpose(JT);
        if(debug){
           console.log("J:");
           print_multi_array(J);
        }
        //J = [
        //    [-0.1,  0.0,  0.0 ], Jvx
        //    [ 0.0, -1.0, -1.0 ], Jvy
        //    [ 0.0,  0.0,  0.0 ], Jvz
        //    [ 0.0,  0.0,  0.0 ], Jwx
        //    [ 0.0,  0.0,  0.0 ], Jwy
        //    [ 1.0,  0.0,  0.0 ]  Jwz
        //]
        
        //Simplify the Jacobian matrix. N<6. Under-actuated robot.
        var J = [
            [J[0][0],J[0][1],J[0][2]], //Jvx
            [J[1][0],J[1][1],J[1][2]], //Jvy
            [J[2][0],J[2][1],J[2][2]]  //Jvz
        ];
        if(debug){
           console.log("Jxyz:");
           print_multi_array(J);
        }
        
        //Inverse. See: centralCamera.js function image_jacobian_inverse().
        // Pseudo inverse (right inverse). J* = J^T * (J*J^T)-1 (MATLAB pinv()).
        //    Part 1. J*J^T (3x3 x (3x3)^T = 3x3 x 3x3 = 3x3)
        var JJT = hlao.matrix_multiplication(J,hlao.matrix_transpose(J));
        //var dim = size(JJT); //3x3
        //var m = dim[0]; //Number of rows.
        //var n = dim[1]; //Number of columns.
        var m; var n;
        [m,n] = [JJT.length,JJT[0].length];
        if(debug) console.log('J x JT: ' + m + ' x ' + n);
        //    Part 2. Adjust the 'JJT' matrix for SVD.
        for(var i=0;i<m;i=i+1){ //3x3 to 3x4. For each row of JJT[] add an extra element '0.0' to the beginning of the array (beginning of the row).
            JJT[i].unshift(0.0);
        }
        var offsetRow = []; //Create an array of zeroes (row vector) 1x4.
        for(var i=0;i<(n+1);i=i+1){ //Only 3 + 1 elements as JJT[] is 3 x 3 matrix or 3 x 4 matrix with padded zeroes.
            offsetRow.push(0.0);
        }
        JJT.unshift(offsetRow); //Add the row vector of zeroes to the beginning of JJT[] array - now a 4 x 4 matrix.
        //    Part 3. SVD.
        var w = hlao.zeros_vector((n+1),'row'); //Row vector where index = 0 is undefined.
        var v = hlao.zeros_matrix((n+1),(n+1)); //Matrix.
        var uwv = svdcmp.svdcmp(JJT, m, n, w, v);
        var U = svdcmp.svdClean(uwv[0]); //Drop the first element in the array as it is zero.
        var S = svdcmp.svdClean(uwv[1]); //W
        var V = svdcmp.svdClean(uwv[2]);
        if(debug){
            console.log('U:');
            print_multi_array(U);
            console.log('S:');
            print_multi_array(S);
            console.log('V:');
            print_multi_array(V);
        }
        //ref:
        //   - http://www.kwon3d.com/theory/jkinem/svd.html
        //   - http://au.mathworks.com/help/matlab/ref/svd.html
        //   - https://au.mathworks.com/help/matlab/ref/pinv.html
        //   - https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse
        var diag = [];
        var TOL = 1e-12; //Tolerance for the evaluation of 'S'.
        for(var i=0;i<3;i=i+1){ //IMPORTANT: Magic number.
            if(Math.abs(S[i]) > (0.0 + TOL)) diag[i] = 1.0/S[i];
            else diag[i] = 0.0;
        }
        var Sinv = [ //IMPORTANT: fix magic number.
            [diag[0],     0.0,     0.0],
            [    0.0, diag[1],     0.0],
            [    0.0,     0.0, diag[2]],
        ];
        var JJTinv = hlao.matrix_multiplication(
                hlao.matrix_multiplication(
                    V,Sinv
                ),
                hlao.matrix_transpose(U)
            );
        if(debug){
            console.log('(J x JT)-1:');
            print_multi_array(JJTinv);
        }
        //Part 4.  J^T * (J*J^T)-1 ((3x3)^T x 3x3 = 3x3)
        var Jstar = hlao.matrix_multiplication(hlao.matrix_transpose(J),JJTinv);
        if(debug){
            console.log('J*:');
            print_multi_array(Jstar);
        }
        //Jstar = [
        //    [-10.0, 0.0,0.0],
        //    [  0.0,-0.5,0.0],
        //    [  0.0,-0.5,0.0],
        //];
        
        //Resolved-Rate Motion Control. Ref: Robotics, Vision and Control, page 180.
        //q_dot_star<k> = J(q<k>)-1 * v_star
        //q_star<i+1> = q<k> + delta_t * q_dot_star<k>
        //where:
        //   - delta_t is the sample interval
        //   - q_star<i+1> the desired joint angles for the next time step
        //   - q<k> the current joint angles
        //   - v_star the desired end-effector velocity
        
        //Joint velocity ('q_dot_star').
        var qd = hlao.matrix_multiplication(Jstar,v_star);
        if(1){
            console.log('qd:');
            print_multi_array(qd);
        }
        //qd = [[-1],[0],[0]];
        
        //Update the joint coordinates ('q_star').
        links[0][3] = links[0][3] + dt*qd[0][0]; //revolute
        links[2][2] = links[2][2] + dt*qd[1][0]; //prismatic
        links[4][2] = links[4][2] + dt*qd[2][0]; //prismatic
        if(debug){
            console.log('links:');
            print_multi_array(links);
        }
    }
    
    //Calculate homogeneous matrix of each link.
    var T = [];
    T[0] = [ //z0
        [1.0,0.0,0.0,0.0],
        [0.0,1.0,0.0,0.0],
        [0.0,0.0,1.0,0.0],
        [0.0,0.0,0.0,1.0]
    ];
    for(var i=0;i<nbr_links;i=i+1){
        T.push(mcht.Aij(links[i])); //transform from from {Ai} to {Ai-1}
        if(debug){
            console.log("T[" + i + "]:");
            print_multi_array(T[i]);
        }
    }
            
    //Calculate homogeneous matrix of each link frame with reference to world frame.
    var T0 = [];
    for(var i=0;i<=nbr_links;i=i+1){
        if(i == 0){
            T0[i] = T[i];
        } else {
            T0[i] = hlao.matrix_multiplication(T0[i-1],T[i]);
        }
        if(debug){
            console.log("T0[" + i + "]:");
            print_multi_array(T0[i]);
        }
    }
    
    //Get the end-effector position.
    var pe = [
        [T0[nbr_links][0][3]],
        [T0[nbr_links][1][3]],
        [T0[nbr_links][2][3]]
    ];
    if(1){
        console.log("pe (final):");
        print_multi_array(pe);
    }
    
    return([qd,T0cam]);
}

function print_multi_array(A){
    var m = A.length;
    var n = A[0].length;
    var str = "";
    
    if(typeof n == 'undefined'){
        //row vector
        m = 1;
        n = A.length;
        str = str + '[';
        for(var j=0;j<n;j=j+1){ //column
            str = str + A[j];
            if(j < n-1) str = str + ',';
        }
        str = str + '];';
    } else {
        str = str + '[\n';
        for(var i=0;i<m;i=i+1){ //row
            str = str + '    [';
            for(var j=0;j<n;j=j+1){ //column
                str = str + A[i][j];
                if(j < n-1) str = str + ',';
            }
            str = str + ']';
            if(i < m-1) str = str + ',';
            str = str + '\n';
        }
        str = str + '];';
    }
    
    console.log(str);
}

export {
    IR_rob3,
    print_multi_array
};
