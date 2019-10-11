function evaluarAltura( altitud ) {
    if( altitud < 500 ) {
        return "hipotesisA";
    } else if ( altitud < 1000 ) {
        return "hipotesisB";
    } else if ( altitud >= 1000 ) {
        return "hipotesisC";
    }

    return "errorLecturaDatos";
}

function cuberoot(x) {
    var y = Math.pow(Math.abs(x), 1/3);
    return x < 0 ? -y : y;
}

function solveCubic(a, b, c, d) {
    if (Math.abs(a) < 1e-8) { // Quadratic case, ax^2+bx+c=0
        a = b; b = c; c = d;
        if (Math.abs(a) < 1e-8) { // Linear case, ax+b=0
            a = b; b = c;
            if (Math.abs(a) < 1e-8) // Degenerate case
                return [];
            return [-b/a];
        }

        var D = b*b - 4*a*c;
        if (Math.abs(D) < 1e-8)
            return [-b/(2*a)];
        else if (D > 0)
            return [(-b+Math.sqrt(D))/(2*a), (-b-Math.sqrt(D))/(2*a)];
        return [];
    }

    // Convert to depressed cubic t^3+pt+q = 0 (subst x = t - b/3a)
    var p = (3*a*c - b*b)/(3*a*a);
    var q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);
    var roots;

    if (Math.abs(p) < 1e-8) { // p = 0 -> t^3 = -q -> t = -q^1/3
        roots = [cuberoot(-q)];
    } else if (Math.abs(q) < 1e-8) { // q = 0 -> t^3 + pt = 0 -> t(t^2+p)=0
        roots = [0].concat(p < 0 ? [Math.sqrt(-p), -Math.sqrt(-p)] : []);
    } else {
        var D = q*q/4 + p*p*p/27;
        if (Math.abs(D) < 1e-8) {       // D = 0 -> two roots
            roots = [-1.5*q/p, 3*q/p];
        } else if (D > 0) {             // Only one real root
            var u = cuberoot(-q/2 - Math.sqrt(D));
            roots = [u - p/(3*u)];
        } else {                        // D < 0, three roots, but needs to use complex numbers/trigonometric solution
            var u = 2*Math.sqrt(-p/3);
            var t = Math.acos(3*q/p/u)/3;  // D < 0 implies p < 0 and acos argument in [-1..1]
            var k = 2*Math.PI/3;
            roots = [u*Math.cos(t), u*Math.cos(t-k), u*Math.cos(t-2*k)];
        }
    }

    // Convert back from depressed cubic
    for (var i = 0; i < roots.length; i++)
        roots[i] -= b/(3*a);

    return roots;
}

function operacionesComuneHipotesisPrincipales( cargaRotura, seccion, pesoConductor, tipoCable ) {
    valores = new Array();

    /*
        Traccion maxima T = CR/3
        Donde:
        CR, es la Carga de rotura
    */
   valores['traccion'] = cargaRotura/( tipoCable == false ? 3:5 );

    /*
        ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        t = T/S
        Donde:
        T, es la Traccion maxima
        S, es la seccion
    */
   valores['t'] = valores.traccion/seccion;

    /*
        Peso por metro y milimetro cuadrado W = Pc/S
        Donde:
        Pc, es el peso del conductor
        S, es la seccion
    */
   valores['W'] = pesoConductor/seccion;

    return valores;
}

function operacionesComuneHipotesisPrincipalesConVariables( vanoPromedio, moduloElasticidad, Pprima, pesoConductor, W, t ) {
    valores = new Array();

    /*
        Coeficiente de sobrecarga inicial M = P'/Pc
        Donde:
        P', ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        Pc, es el peso del conductor
    */
   valores['M'] =Pprima/pesoConductor;

    /*
        Constante K en la ecuacion de cambio de condiciones k = t - [ ( a^2*M^2*W^2*E )/( 24*t^2 ) ]
        Donde:
        t, ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        a, es el Vano promedio
        M, es el coeficiente de sobrecarga inicial
        W, es el peso por metro y milimetro cuadrado
        E, es el modulo de elasticidad
    */
   valores['k'] = t - ( ( Math.pow(vanoPromedio, 2)*Math.pow(valores.M, 2)*Math.pow(W, 2)*moduloElasticidad )/( 24*Math.pow(t, 2) ) );

    /*
        Flecha Vertical (No hay viento) Fc = ( a^2*W*M )/ (8*t)
        Donde:
        t, ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        a, es el Vano promedio
        M, es el coeficiente de sobrecarga inicial
        W, es el peso por metro y milimetro cuadrado
    */
    valores['Fc'] = ( Math.pow(vanoPromedio, 2)*W*valores.M )/( 8*t );

   return valores;
}

function calculoSobrecarga( pesoConductor, diametro ) {
    var P = diametro <= 16 ? 60*diametro/1000:50*diametro/1000;
    var Pprima = Math.sqrt( Math.pow(pesoConductor, 2) + Math.pow(P, 2) );
    
    var M = Pprima/pesoConductor; 

    return M;
}

function calculoTension( diferencialTemperatura, k, coeficienteDilatacion, moduloElasticidad, vanoPromedio, W, sobrecarga ) {
    var tension = null;

    var x = k - coeficienteDilatacion*moduloElasticidad*diferencialTemperatura;
    var y = ( Math.pow(vanoPromedio, 2)*Math.pow(sobrecarga, 2)*Math.pow(W, 2)*moduloElasticidad )/24;
    
    // Resuelve la ecuacion cubica y busca el valor real
    solveCubic( 1, -x, 0 , -y ).forEach(element => {
        if(element > 0) { tension = element}
    });

    return tension;
}

// Hipotesis
function hipotesisA( cargaRotura, seccion, diametro, pesoConductor, vanoPromedio, moduloElasticidad, tipoCable ) {
    var respuesta = new Array();
    respuesta = operacionesComuneHipotesisPrincipales( cargaRotura, seccion, pesoConductor, tipoCable );
    
    respuesta['temperatura'] = -5;

    /*
        ##FATA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        P = 0.50*d/1000
        Donde:
        d, es el Diametro
    */
    respuesta['P'] = diametro <= 16 ? 60*diametro/1000:50*diametro/1000;

    /*
        ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        P' = sqrt ( P^2 + Pc^2 )
        Donde:
        P, es ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        Pc, es el peso del conductor
    */
    respuesta['Pprima'] = Math.sqrt( Math.pow(respuesta.P, 2) + Math.pow(pesoConductor, 2) );
    
    valores = operacionesComuneHipotesisPrincipalesConVariables( vanoPromedio, moduloElasticidad, respuesta.Pprima, pesoConductor, respuesta.W, respuesta.t );
    
    for ( valor in valores ) {
        respuesta[valor] = valores[valor];
    }

    return respuesta;
}

function hipotesisB( cargaRotura, seccion, diametro, pesoConductor, vanoPromedio, moduloElasticidad, tipoCable ) {
    var respuesta = new Array();
    respuesta = operacionesComuneHipotesisPrincipales( cargaRotura, seccion, pesoConductor, tipoCable );

    respuesta['temperatura'] = -15;

    /*
        ##FATA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        P = 0.18*sqrt(d)
        Donde:
        d, es el Diametro
    */
    respuesta['P'] = 0.18*Math.sqrt(diametro);

    /*
        ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        P' = P + Pc
        Donde:
        P, es ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        Pc, es el peso del conductor
    */
    respuesta['Pprima'] = respuesta.P + pesoConductor;
    
    valores = operacionesComuneHipotesisPrincipalesConVariables( vanoPromedio, moduloElasticidad, respuesta.Pprima, pesoConductor, respuesta.W, respuesta.t );
    
    for ( valor in valores ) {
        respuesta[valor] = valores[valor];
    }

    return respuesta;
}

function hipotesisC( cargaRotura, seccion, diametro, pesoConductor, vanoPromedio, moduloElasticidad, tipoCable ) {
    var respuesta = new Array();
    respuesta = operacionesComuneHipotesisPrincipales( cargaRotura, seccion, pesoConductor, tipoCable );

    respuesta['temperatura'] = -20;

    /*
        ##FATA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        P = 0.36*sqrt(d)
        Donde:
        d, es el Diametro
    */
    respuesta['P'] = 0.36*Math.sqrt(diametro);

    /*
        ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        P' = P + Pc
        Donde:
        P, es ##FALTA LA DESCRIPCION DE COMO SE LLAMA ESTE VALOR
        Pc, es el peso del conductor
    */
    respuesta['Pprima'] = respuesta.P + pesoConductor;
    
    valores = operacionesComuneHipotesisPrincipalesConVariables( vanoPromedio, moduloElasticidad, respuesta.Pprima, pesoConductor, respuesta.W, respuesta.t );
    
    for ( valor in valores ) {
        respuesta[valor] = valores[valor];
    }

    return respuesta;
}

function hipotesisSecuenciales( temperaturaX, temperatura, sobrecarga, diametro, pesoConductor, k, coeficienteDilatacion, moduloElasticidad, vanoPromedio, W, seccion, cargaRotura ) {
    var respuesta = new Array();
    var diferencialTemperatura = temperatura - temperaturaX;

    respuesta['temperatura'] = temperatura;
    respuesta['diferencialTemperatura'] = diferencialTemperatura;
    respuesta['M'] = sobrecarga;
    respuesta['t'] = calculoTension( diferencialTemperatura, k, coeficienteDilatacion, moduloElasticidad, vanoPromedio, W, sobrecarga );
    respuesta['traccion'] = respuesta.t*seccion;
    respuesta['coeficienteSeguridad'] = cargaRotura/respuesta.traccion;
    respuesta['flechaInclinada'] = ( Math.pow(vanoPromedio, 2)*W*sobrecarga )/( 8*respuesta.t );

    return respuesta;
}

function llamarHipotesis( cargaRotura, seccion, diametro, pesoConductor, vanoPromedio, moduloElasticidad, altitud, coeficienteDilatacion, tipoCable ) {
    // Arreglo que contendra la informacion resultante
    var datos = new Array();

    // Hipotesis Primaria
    primeraHipotesis = evaluarAltura( altitud );
    datos['X'] = window[primeraHipotesis](
        cargaRotura,
        seccion,
        diametro,
        pesoConductor,
        vanoPromedio,
        moduloElasticidad,
        tipoCable
    );

    datos.hipotesis = primeraHipotesis.substring(9);
    
    // Hipotesis secuenciales
    // Loop para ir  de D a J
	for(i = 68; i < 75; i++){
		// Convertir el codigo CHAR a STRING
        var letra =String.fromCharCode(i);
        var temperatura = null;
        var sobrecarga = 1;

        if( letra === 'D' || letra === 'G' ) {
            temperatura = 15;

            if( letra === 'D') {
                sobrecarga = calculoSobrecarga(pesoConductor, diametro);
            }
        } else if( letra === 'E' ) {
            temperatura = 50;
        } else if( letra === 'F' ) {
            temperatura = 0;

            sobrecarga = (primeraHipotesis === 'hipotesisA' ? 1:datos.X.M);
        } else if( letra === 'G' ) {
            temperatura = 15;
        } else if( letra === 'H' || letra === 'J' ) {
            temperatura = -5;

            if( letra === 'J') {
                sobrecarga = calculoSobrecarga(pesoConductor, diametro);
            }
        } else if( letra === 'I' ) {
            temperatura = datos.X.temperatura;
        }

        datos[letra] = hipotesisSecuenciales(
            datos.X.temperatura,
            temperatura,
            sobrecarga,
            diametro,
            pesoConductor,
            datos.X.k,
            coeficienteDilatacion,
            moduloElasticidad,
            vanoPromedio,
            datos.X.W,
            seccion,
            cargaRotura
        );
    }
    
    return datos;
}

function main() {
    // Caracteristicas del conductor
    var diametro = parseFloat(document.getElementById('diametro').value,4); // Abreviacion en documentos 'd'
    var seccion = parseFloat(document.getElementById('seccion').value,4); // Abreviacion en documentos 'S'
    var pesoConductor = parseFloat(document.getElementById('pesoConductor').value,4); // Abreviacion en documentos 'Pc'
    var cargaRotura = parseFloat(document.getElementById('cargaRotura').value,4); // Abreviacion en documentos 'CR'
    var moduloElasticidad = parseFloat(document.getElementById('moduloElasticidad').value,4); // Abreviacion en documentos 'E'
    var coeficienteDilatacion = parseFloat(document.getElementById('coeficienteDilatacion').value,10); // Abreviacion en documentos 'alfa'

    // Caracteristicas de la linea
    var altitud = parseFloat(document.getElementById('altitud').value,4); // Abreviacion en documentos 'h'
    var vanoPromedio = parseFloat(document.getElementById('vano').value,4); // Abreviacion en documentos 'a'
    var longitud = null; // Abreviacion en documentos 'l'
    var voltaje = null; // Abreviacion en documentos 'U'

    var tipoCable = document.querySelector('input[name="tipoCable"]:checked').value == 'tierra' ? true:false;

    var datos = llamarHipotesis(
        cargaRotura,
        seccion,
        diametro,
        pesoConductor,
        vanoPromedio,
        moduloElasticidad,
        altitud,
        coeficienteDilatacion,
        tipoCable
    );

    console.log(datos);

    //HIPOTESIS X
    document.getElementById('hipotesisX').innerHTML = datos.hipotesis;
    document.getElementById('temperaturaHipotesisX').innerHTML = datos.X.temperatura + ( datos.hipotesis == 'A' ? '°C. Viento':'°C. Hielo' );
    document.getElementById('tensionX').innerHTML = parseFloat(datos.X.traccion).toFixed(4);
    document.getElementById('coeficienteX').innerHTML = ( tipoCable == false ? 3:5 );
    document.getElementById('flechaVX').innerHTML =  ( datos.hipotesis == 'A' ? '-':parseFloat(datos.X.Fc).toFixed(4) );
    document.getElementById('flechaIX').innerHTML = ( datos.hipotesis == 'A' ? parseFloat(datos.X.Fc).toFixed(4):'-' );

    //HIPOTESIS D
    document.getElementById('tensionD').innerHTML = parseFloat(datos.D.traccion).toFixed(2);
    document.getElementById('coeficienteD').innerHTML = parseFloat(datos.D.coeficienteSeguridad).toFixed(4);
    document.getElementById('flechaVD').innerHTML = '-';
    document.getElementById('flechaID').innerHTML = parseFloat(datos.D.flechaInclinada).toFixed(4);

    //HIPOTESIS E
    document.getElementById('tensionE').innerHTML = parseFloat(datos.E.traccion).toFixed(2);
    document.getElementById('coeficienteE').innerHTML = parseFloat(datos.E.coeficienteSeguridad).toFixed(4);
    document.getElementById('flechaVE').innerHTML = parseFloat(datos.E.flechaInclinada).toFixed(4);
    document.getElementById('flechaIE').innerHTML = '-';

    //HIPOTESIS F
    document.getElementById('tensionF').innerHTML = parseFloat(datos.F.traccion).toFixed(2);
    document.getElementById('coeficienteF').innerHTML = parseFloat(datos.F.coeficienteSeguridad).toFixed(4);
    document.getElementById('flechaVF').innerHTML = parseFloat(datos.F.flechaInclinada).toFixed(4);
    document.getElementById('flechaIF').innerHTML = '-';
    
    //HIPOTESIS G
    document.getElementById('tensionG').innerHTML = parseFloat(datos.G.traccion).toFixed(2);
    document.getElementById('coeficienteG').innerHTML = parseFloat(datos.G.coeficienteSeguridad).toFixed(4);
    document.getElementById('flechaVG').innerHTML = parseFloat(datos.G.flechaInclinada).toFixed(4);
    document.getElementById('flechaIG').innerHTML = '-';

    //HIPOTESIS H
    document.getElementById('tensionH').innerHTML = parseFloat(datos.H.traccion).toFixed(4);
    document.getElementById('coeficienteH').innerHTML = parseFloat(datos.H.coeficienteSeguridad).toFixed(4);
    document.getElementById('flechaVH').innerHTML = parseFloat(datos.H.flechaInclinada).toFixed(4);
    document.getElementById('flechaIH').innerHTML = '-';

    //HIPOTESIS I
    document.getElementById('tensionI').innerHTML = parseFloat(datos.I.traccion).toFixed(4);
    document.getElementById('coeficienteI').innerHTML = parseFloat(datos.I.coeficienteSeguridad).toFixed(4);
    document.getElementById('flechaVI').innerHTML = parseFloat(datos.I.flechaInclinada).toFixed(4);
    document.getElementById('flechaII').innerHTML = '-';

    //HIPOTESIS J
    document.getElementById('tensionJ').innerHTML = parseFloat(datos.J.traccion).toFixed(4);
    document.getElementById('coeficienteJ').innerHTML = parseFloat(datos.J.coeficienteSeguridad).toFixed(4);
    document.getElementById('flechaVJ').innerHTML = '-';
    document.getElementById('flechaIJ').innerHTML = parseFloat(datos.J.flechaInclinada).toFixed(4);

    //TCD - THF
    document.getElementById('tcd').innerHTML = parseFloat((datos.G.traccion/cargaRotura)*100).toFixed(2) + '%';
    document.getElementById('thf').innerHTML = parseFloat((datos.H.traccion/cargaRotura)*100).toFixed(2) + '%';

    //OBSERVACIONES
    //Flecha maxima vertical
    var maximaVertical = Math.max(
        ( datos.hipotesis == 'A' ? null:datos.X.Fc ),
        datos.E.flechaInclinada,
        datos.F.flechaInclinada,
        datos.G.flechaInclinada,
        datos.H.flechaInclinada,
        datos.I.flechaInclinada
    );

    //Flecha maxima inclinada
    var maximaInclinada = Math.max(
        ( datos.hipotesis == 'A' ? datos.X.Fc:null ),
        datos.D.flechaInclinada,
        datos.J.flechaInclinada
    );

    //Flecha minima
    var minimaVertical = Math.min(
        ( datos.hipotesis == 'A' ? null:datos.X.Fc ),
        datos.E.flechaInclinada,
        datos.F.flechaInclinada,
        datos.G.flechaInclinada,
        datos.H.flechaInclinada,
        datos.I.flechaInclinada
    );

    //Coeficiente de seguridad minimo
    var minimaCoeficienteSeguridad = Math.min(
        3,
        datos.D.coeficienteSeguridad,
        datos.E.coeficienteSeguridad,
        datos.F.coeficienteSeguridad,
        datos.G.coeficienteSeguridad,
        datos.H.coeficienteSeguridad,
        datos.I.coeficienteSeguridad,
        datos.J.coeficienteSeguridad
    );

    //Coeficiente de seguridad maximo
    var maximoCoeficienteSeguridad = Math.max(
        3,
        datos.D.coeficienteSeguridad,
        datos.E.coeficienteSeguridad,
        datos.F.coeficienteSeguridad,
        datos.G.coeficienteSeguridad,
        datos.H.coeficienteSeguridad,
        datos.I.coeficienteSeguridad,
        datos.J.coeficienteSeguridad
    );

    document.getElementById('observacionX').innerHTML = "<p> Tracci&oacute;n m&aacute;xima"+ ( ( maximaInclinada && datos.hipotesis == 'A' ) == datos.D.flechaInclinada ? 'Flecha m&aacute;xima inclinada<br>':'') + ( ( maximaVertical && datos.hipotesis != 'A' ) == datos.X.Fc ? 'Flecha m&aacute;xima vertical<br>':'') + (minimaVertical == datos.X.Fc ? 'Flecha m&iacute;nima vertical<br>':'') + (maximoCoeficienteSeguridad == 3 ? 'Coeficiente de seguridad m&aacute;ximo<br>':'') + (minimaCoeficienteSeguridad == 3 ? '<br>Coeficiente de seguridad m&iacute;nimo':'') + "</p> ";
    document.getElementById('observacionD').innerHTML = "<p>"+ (maximaInclinada == datos.D.flechaInclinada ? 'Flecha m&aacute;xima inclinada<br>':'') + (maximoCoeficienteSeguridad == datos.D.coeficienteSeguridad ? 'Coeficiente de seguridad m&aacute;ximo<br>':'') + (minimaCoeficienteSeguridad == datos.D.coeficienteSeguridad ? 'Coeficiente de seguridad m&iacute;nimo':'') + "</p> ";
    document.getElementById('observacionE').innerHTML = "<p>"+ (maximaVertical == datos.E.flechaInclinada ? 'Flecha m&aacute;xima vertical<br>':'') + (minimaVertical == datos.E.flechaInclinada ? 'Flecha m&iacute;nima vertical<br>':'') + (maximoCoeficienteSeguridad == datos.E.coeficienteSeguridad ? 'Coeficiente de seguridad m&aacute;ximo<br>':'') + (minimaCoeficienteSeguridad == datos.E.coeficienteSeguridad ? '<br>Coeficiente de seguridad m&iacute;nimo':'') + "</p> ";
    document.getElementById('observacionF').innerHTML = "<p>"+ (maximaVertical == datos.F.flechaInclinada ? 'Flecha m&aacute;xima vertical<br>':'') + (minimaVertical == datos.F.flechaInclinada ? 'Flecha m&iacute;nima vertical<br>':'') + (maximoCoeficienteSeguridad == datos.F.coeficienteSeguridad ? 'Coeficiente de seguridad m&aacute;ximo<br>':'') + (minimaCoeficienteSeguridad == datos.F.coeficienteSeguridad ? '<br>Coeficiente de seguridad m&iacute;nimo':'') + "</p> ";
    document.getElementById('observacionG').innerHTML = "<p>"+ (maximaVertical == datos.G.flechaInclinada ? 'Flecha m&aacute;xima vertical<br>':'') + (minimaVertical == datos.G.flechaInclinada ? 'Flecha m&iacute;nima vertical<br>':'') + (maximoCoeficienteSeguridad == datos.G.coeficienteSeguridad ? 'Coeficiente de seguridad m&aacute;ximo<br>':'') + (minimaCoeficienteSeguridad == datos.G.coeficienteSeguridad ? '<br>Coeficiente de seguridad m&iacute;nimo':'') + "</p> ";
    document.getElementById('observacionH').innerHTML = "<p>"+ (maximaVertical == datos.H.flechaInclinada ? 'Flecha m&aacute;xima vertical<br>':'') + (minimaVertical == datos.H.flechaInclinada ? 'Flecha m&iacute;nima vertical<br>':'') + (maximoCoeficienteSeguridad == datos.H.coeficienteSeguridad ? 'Coeficiente de seguridad m&aacute;ximo<br>':'') + (minimaCoeficienteSeguridad == datos.H.coeficienteSeguridad ? '<br>Coeficiente de seguridad m&iacute;nimo':'') + "</p> ";
    document.getElementById('observacionI').innerHTML = "<p>"+ (maximaVertical == datos.I.flechaInclinada ? 'Flecha m&aacute;xima vertical<br>':'') + (minimaVertical == datos.I.flechaInclinada ? 'Flecha m&iacute;nima vertical<br>':'') + (maximoCoeficienteSeguridad == datos.I.coeficienteSeguridad ? 'Coeficiente de seguridad m&aacute;ximo<br>':'') + (minimaCoeficienteSeguridad == datos.I.coeficienteSeguridad ? '<br>Coeficiente de seguridad m&iacute;nimo':'') + "</p> ";
    document.getElementById('observacionJ').innerHTML = "<p>"+ (maximaInclinada == datos.J.flechaInclinada ? 'Flecha m&aacute;xima inclinada<br>':'') + (maximoCoeficienteSeguridad == datos.J.coeficienteSeguridad ? 'Coeficiente de seguridad m&aacute;ximo<br>':'') + (minimaCoeficienteSeguridad == datos.J.coeficienteSeguridad ? 'Coeficiente de seguridad m&iacute;nimo':'') + "</p> ";
}

function limpiar() {
    window.location.reload();
}