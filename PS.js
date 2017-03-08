<!--Define JavaScript functions.-->
// this code in part adapted from f1z.js, taken from http://www.akiti.ca/
// by Tony Withers, 5 Sept 2004.
var PScoeff = new Array(10);
for (var i = 0; i < 10; i++) {
	PScoeff[i] = new Array(6);
	for (var j = 0; j < 10; j++) PScoeff[i][j] =0;
}

PScoeff[0][2]=0.24657688*Math.pow(10,6) ;
PScoeff[0][3]=0.51359951*Math.pow(10,2) ;
PScoeff[1][2]=0.58638965*Math.pow(10,0) ;
PScoeff[1][3]=-0.28646939*Math.pow(10,-2) ;
PScoeff[1][4]=0.31375577*Math.pow(10,-4) ;
PScoeff[2][2]=-0.62783840*Math.pow(10,1) ;
PScoeff[2][3]=0.14791599*Math.pow(10,-1) ;
PScoeff[2][4]=0.35779579*Math.pow(10,-3) ;
PScoeff[2][5]=0.15432925*Math.pow(10,-7) ;
PScoeff[3][3]=-0.42719875*Math.pow(10,0) ;
PScoeff[3][4]=-0.16325155*Math.pow(10,-4) ;
PScoeff[4][2]=0.56654978*Math.pow(10,4) ;
PScoeff[4][3]=-0.16580167*Math.pow(10,2) ;
PScoeff[4][4]=0.76560762*Math.pow(10,-1) ;
PScoeff[5][3]=0.10917883*Math.pow(10,0) ;
PScoeff[6][0]=0.38878656*Math.pow(10,13) ;
PScoeff[6][1]=-0.13494878*Math.pow(10,9) ;
PScoeff[6][2]=0.30916564*Math.pow(10,6) ;
PScoeff[6][3]=0.75591105*Math.pow(10,1) ;
PScoeff[7][2]=-0.65537898*Math.pow(10,5) ;
PScoeff[7][3]=0.18810675*Math.pow(10,3) ;
PScoeff[8][0]=-0.14182435*Math.pow(10,14) ;
PScoeff[8][1]=0.18165390*Math.pow(10,9) ;
PScoeff[8][2]=-0.19769068*Math.pow(10,6) ;
PScoeff[8][3]=-0.23530318*Math.pow(10,2) ;
PScoeff[9][2]=0.92093375*Math.pow(10,5) ;
PScoeff[9][3]=0.12246777*Math.pow(10,3) ;
var cs= new Array(10);

function eos(T,V) {
	// takes T(K) and V(cc/mol)	
	var den=1/V;
	var R=8314472;
	
	var var_num, var_denom;
	
	var_num=cs[2]+2*cs[3]*den+3*cs[4]*Math.pow(den,2)+4*cs[5]*Math.pow(den,3);
	var_denom=Math.pow((cs[1]+cs[2]*den+cs[3]*Math.pow(den,2)+cs[4]*Math.pow(den,3)+cs[5]*Math.pow(den,4)),2);
	pressure=den+cs[0]*Math.pow(den,2)-Math.pow(den,2)*(var_num/var_denom);
	pressure+=cs[6]*Math.pow(den,2)*Math.exp(-cs[7]*den)+cs[8]*Math.pow(den,2)*Math.exp(-cs[9]*den);
	pressure*=R*T;
	return pressure;		// pressure in Pa
}

function PSfug(P,T,V) {
	var den=1/V;
	var lnf, quotient;
	var R=8314472;
	
	quotient = cs[0]*den+(1/(cs[1]+cs[2]*den+cs[3]*Math.pow(den,2)+cs[4]*Math.pow(den,3)+cs[5]*Math.pow(den,4))-1/cs[1]);
	quotient-= cs[6]/cs[7]*(Math.exp(-cs[7]*den)-1);
	quotient-= cs[8]/cs[9]*(Math.exp(-cs[9]*den)-1);
	lnf=(Math.log(den)+quotient+P/(den*R*T))+Math.log(R*T)-1;
	return Math.exp(lnf)/1e9; // fugacity in GPa	
}

function sign(y){
 return ((y < 0) ? -1 : 1);
}

function findroot(form){
var pressure = parseFloat(form.pressure.value)*1000000000; //Pa
var temperature = parseFloat(form.temperature.value)+273.15; // K
var b=parseFloat(form.lowvol.value); // cc/mol
var c=parseFloat(form.highvol.value); // cc/mol

for (i = 0; i < 10; i++) cs[i]=PScoeff[i][0]*Math.pow(temperature,-4)+PScoeff[i][1]*Math.pow(temperature,-2)+PScoeff[i][2]*Math.pow(temperature,-1)+PScoeff[i][3]+PScoeff[i][4]*temperature+PScoeff[i][5]*Math.pow(temperature,2);
	
if (pressure <= 0)  {
	alert("Pressure out of range");
	return;
}

var z = (b + c)/2;
var fz = eos(temperature,z)-pressure;
var fc = fz, t = b, count = 2;
var fb = eos(temperature,b)-pressure;

if (fb == 0){
 form.volume.value = b;
 //alert("A root has been found, but the interval has not collapsed to the requested tolerance");
 form.fugacity.value=PSfug(pressure,temperature,b);
 //form.funcEval.value = count;
 return;
}

if (fz == 0){
 form.volume.value = z;
 //alert("A root has been found, but the interval has not collapsed to the requested tolerance");
 form.fugacity.value=PSfug(pressure,temperature,z);
 //form.funcEval.value = count;
 return;
}

if (sign(fz) == sign(fb))
{t = c;
 fc = eos(temperature,c)-pressure;
 count = 3;

 if (fc == 0){
  form.volume.value = c;
  //alert("A root has been found, but the interval has not collapsed to the requested tolerance");
  form.fugacity.value=PSfug(pressure,temperature,c);
  //form.funcEval.value = count;
  return;
 }

 if (sign(fz) != sign(fc))
 {b = z;
  fb = fz;
 }
 else {//Sign is the same on this interval as well; no zero on input interval
  alert("volume not bracketed");
  return;
 }
}
else c = z;

var a = c, fa = fc, acmb, ic = 0, tol, p, q, acbs = Math.abs(-b + c), maxit = 1000, ae, cmb;

acmb = 1.0;
do {
  ae = acmb;
  acmb /= 2.0;
  cmb = 1.0 + acmb;
} while (cmb > 1.0);
// ae should equal the machine epsilon
// form.epmch.value = ae;

do {
 if (Math.abs(fc) < Math.abs(fb))  //Interchange if necessary.
 {a = b;
  fa = fb;
  b = c;
  fb = fc;
  c = a;
  fc = fa
 }
 cmb = (-b + c)/2;
 acmb = Math.abs(cmb);
 tol = Math.abs(b);
 tol++;
 tol *= ae;

 if (acmb <= tol) break;

 if (fb == 0) 
 {//alert("A root has been found, but the interval has not collapsed to the requested tolerance");
  break;
 }

// Calculate new iterate implicitly as b + p/q, where p is arranged to be >= 0. 
// This implicit form is used to prevent overflow.

 p = (-a + b)*fb;
 q = -fb + fa;
 if (p < 0)
 {p = -p;
  q = -q;
 }

// Update a and check for satisfactory reduction in the size of the bracketing interval. 
// If not, perform bisection.

 a = b;
 fa = fb;
 ic++;

 if ((ic >= 4) && (8*acmb >= acbs))
  b = (c + b)/2;  //Use bisection
 else
 {if (ic >= 4)
  {ic = 0;
   acbs = acmb;
  }
  if (p <= tol*Math.abs(q))  //Test for too small a change
   b += tol*sign(cmb);
  else    //Root between b and (b + c)/2
  {if (p < cmb*q)  //Use secant rule
    b += p/q;
   else            //Use bisection
    b = (c + b)/2;
  }
 }

// Have now computed new iterate, b.

 fb = eos(temperature,b)-pressure;
 count++;

 if (fb == 0){
  form.volume.value = b;
  //alert("A root has been found, but the interval has not collapsed to the requested tolerance");
  form.fugacity.value=PSfug(pressure,temperature,b);
  //form.funcEval.value = count;
  return;
 }

//Decide whether next step interpolation or extrapolation

 if (sign(fb) == sign(fc))
 {c = a;
  fc = fa;
 } 

} while (count < maxit); //End do-while loop
if (count >= maxit) alert("maximum number of iterations exceeded (more than "+maxit+" iterations)");
form.volume.value = b;
form.fugacity.value=PSfug(pressure,temperature,b);
//form.funcEval.value = count;

return;
}

// end of JavaScript function definitions -->