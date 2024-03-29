var app = angular.module('farmdoc', ['ui.bootstrap', 'ui.router', 'ngCookies', 'ngResource']);
app.constant('config', {
    serverRoot: FARMDOC_CONFIG.serviceUrl
});
app.factory('nav', [function() {
    var _service = {
        _pageName: ""
    };
    _service.setPageName = function(pageName) {
        this._pageName = pageName;
    };
    _service.getPageName = function() {
        return this._pageName;
    }
    return _service;
}]);
"use strict";
var numeric = typeof exports === "undefined" ? function numeric() {} : exports;
if (typeof global !== "undefined") {
    global.numeric = numeric
}
numeric.version = "1.2.6";
numeric.bench = function bench(f, interval) {
    var t1, t2, n, i;
    if (typeof interval === "undefined") {
        interval = 15
    }
    n = .5;
    t1 = new Date;
    while (1) {
        n *= 2;
        for (i = n; i > 3; i -= 4) {
            f();
            f();
            f();
            f()
        }
        while (i > 0) {
            f();
            i--
        }
        t2 = new Date;
        if (t2 - t1 > interval) break
    }
    for (i = n; i > 3; i -= 4) {
        f();
        f();
        f();
        f()
    }
    while (i > 0) {
        f();
        i--
    }
    t2 = new Date;
    return 1e3 * (3 * n - 1) / (t2 - t1)
};
numeric._myIndexOf = function _myIndexOf(w) {
    var n = this.length,
        k;
    for (k = 0; k < n; ++k)
        if (this[k] === w) return k;
    return -1
};
numeric.myIndexOf = Array.prototype.indexOf ? Array.prototype.indexOf : numeric._myIndexOf;
numeric.Function = Function;
numeric.precision = 4;
numeric.largeArray = 50;
numeric.prettyPrint = function prettyPrint(x) {
    function fmtnum(x) {
        if (x === 0) {
            return "0"
        }
        if (isNaN(x)) {
            return "NaN"
        }
        if (x < 0) {
            return "-" + fmtnum(-x)
        }
        if (isFinite(x)) {
            var scale = Math.floor(Math.log(x) / Math.log(10));
            var normalized = x / Math.pow(10, scale);
            var basic = normalized.toPrecision(numeric.precision);
            if (parseFloat(basic) === 10) {
                scale++;
                normalized = 1;
                basic = normalized.toPrecision(numeric.precision)
            }
            return parseFloat(basic).toString() + "e" + scale.toString()
        }
        return "Infinity"
    }
    var ret = [];

    function foo(x) {
        var k;
        if (typeof x === "undefined") {
            ret.push(Array(numeric.precision + 8).join(" "));
            return false
        }
        if (typeof x === "string") {
            ret.push('"' + x + '"');
            return false
        }
        if (typeof x === "boolean") {
            ret.push(x.toString());
            return false
        }
        if (typeof x === "number") {
            var a = fmtnum(x);
            var b = x.toPrecision(numeric.precision);
            var c = parseFloat(x.toString()).toString();
            var d = [a, b, c, parseFloat(b).toString(), parseFloat(c).toString()];
            for (k = 1; k < d.length; k++) {
                if (d[k].length < a.length) a = d[k]
            }
            ret.push(Array(numeric.precision + 8 - a.length).join(" ") + a);
            return false
        }
        if (x === null) {
            ret.push("null");
            return false
        }
        if (typeof x === "function") {
            ret.push(x.toString());
            var flag = false;
            for (k in x) {
                if (x.hasOwnProperty(k)) {
                    if (flag) ret.push(",\n");
                    else ret.push("\n{");
                    flag = true;
                    ret.push(k);
                    ret.push(": \n");
                    foo(x[k])
                }
            }
            if (flag) ret.push("}\n");
            return true
        }
        if (x instanceof Array) {
            if (x.length > numeric.largeArray) {
                ret.push("...Large Array...");
                return true
            }
            var flag = false;
            ret.push("[");
            for (k = 0; k < x.length; k++) {
                if (k > 0) {
                    ret.push(",");
                    if (flag) ret.push("\n ")
                }
                flag = foo(x[k])
            }
            ret.push("]");
            return true
        }
        ret.push("{");
        var flag = false;
        for (k in x) {
            if (x.hasOwnProperty(k)) {
                if (flag) ret.push(",\n");
                flag = true;
                ret.push(k);
                ret.push(": \n");
                foo(x[k])
            }
        }
        ret.push("}");
        return true
    }
    foo(x);
    return ret.join("")
};
numeric.parseDate = function parseDate(d) {
    function foo(d) {
        if (typeof d === "string") {
            return Date.parse(d.replace(/-/g, "/"))
        }
        if (!(d instanceof Array)) {
            throw new Error("parseDate: parameter must be arrays of strings")
        }
        var ret = [],
            k;
        for (k = 0; k < d.length; k++) {
            ret[k] = foo(d[k])
        }
        return ret
    }
    return foo(d)
};
numeric.parseFloat = function parseFloat_(d) {
    function foo(d) {
        if (typeof d === "string") {
            return parseFloat(d)
        }
        if (!(d instanceof Array)) {
            throw new Error("parseFloat: parameter must be arrays of strings")
        }
        var ret = [],
            k;
        for (k = 0; k < d.length; k++) {
            ret[k] = foo(d[k])
        }
        return ret
    }
    return foo(d)
};
numeric.parseCSV = function parseCSV(t) {
    var foo = t.split("\n");
    var j, k;
    var ret = [];
    var pat = /(([^'",]*)|('[^']*')|("[^"]*")),/g;
    var patnum = /^\s*(([+-]?[0-9]+(\.[0-9]*)?(e[+-]?[0-9]+)?)|([+-]?[0-9]*(\.[0-9]+)?(e[+-]?[0-9]+)?))\s*$/;
    var stripper = function(n) {
        return n.substr(0, n.length - 1)
    };
    var count = 0;
    for (k = 0; k < foo.length; k++) {
        var bar = (foo[k] + ",").match(pat),
            baz;
        if (bar.length > 0) {
            ret[count] = [];
            for (j = 0; j < bar.length; j++) {
                baz = stripper(bar[j]);
                if (patnum.test(baz)) {
                    ret[count][j] = parseFloat(baz)
                } else ret[count][j] = baz
            }
            count++
        }
    }
    return ret
};
numeric.toCSV = function toCSV(A) {
    var s = numeric.dim(A);
    var i, j, m, n, row, ret;
    m = s[0];
    n = s[1];
    ret = [];
    for (i = 0; i < m; i++) {
        row = [];
        for (j = 0; j < m; j++) {
            row[j] = A[i][j].toString()
        }
        ret[i] = row.join(", ")
    }
    return ret.join("\n") + "\n"
};
numeric.getURL = function getURL(url) {
    var client = new XMLHttpRequest;
    client.open("GET", url, false);
    client.send();
    return client
};
numeric.imageURL = function imageURL(img) {
    function base64(A) {
        var n = A.length,
            i, x, y, z, p, q, r, s;
        var key = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=";
        var ret = "";
        for (i = 0; i < n; i += 3) {
            x = A[i];
            y = A[i + 1];
            z = A[i + 2];
            p = x >> 2;
            q = ((x & 3) << 4) + (y >> 4);
            r = ((y & 15) << 2) + (z >> 6);
            s = z & 63;
            if (i + 1 >= n) {
                r = s = 64
            } else if (i + 2 >= n) {
                s = 64
            }
            ret += key.charAt(p) + key.charAt(q) + key.charAt(r) + key.charAt(s)
        }
        return ret
    }

    function crc32Array(a, from, to) {
        if (typeof from === "undefined") {
            from = 0
        }
        if (typeof to === "undefined") {
            to = a.length
        }
        var table = [0, 1996959894, 3993919788, 2567524794, 124634137, 1886057615, 3915621685, 2657392035, 249268274, 2044508324, 3772115230, 2547177864, 162941995, 2125561021, 3887607047, 2428444049, 498536548, 1789927666, 4089016648, 2227061214, 450548861, 1843258603, 4107580753, 2211677639, 325883990, 1684777152, 4251122042, 2321926636, 335633487, 1661365465, 4195302755, 2366115317, 997073096, 1281953886, 3579855332, 2724688242, 1006888145, 1258607687, 3524101629, 2768942443, 901097722, 1119000684, 3686517206, 2898065728, 853044451, 1172266101, 3705015759, 2882616665, 651767980, 1373503546, 3369554304, 3218104598, 565507253, 1454621731, 3485111705, 3099436303, 671266974, 1594198024, 3322730930, 2970347812, 795835527, 1483230225, 3244367275, 3060149565, 1994146192, 31158534, 2563907772, 4023717930, 1907459465, 112637215, 2680153253, 3904427059, 2013776290, 251722036, 2517215374, 3775830040, 2137656763, 141376813, 2439277719, 3865271297, 1802195444, 476864866, 2238001368, 4066508878, 1812370925, 453092731, 2181625025, 4111451223, 1706088902, 314042704, 2344532202, 4240017532, 1658658271, 366619977, 2362670323, 4224994405, 1303535960, 984961486, 2747007092, 3569037538, 1256170817, 1037604311, 2765210733, 3554079995, 1131014506, 879679996, 2909243462, 3663771856, 1141124467, 855842277, 2852801631, 3708648649, 1342533948, 654459306, 3188396048, 3373015174, 1466479909, 544179635, 3110523913, 3462522015, 1591671054, 702138776, 2966460450, 3352799412, 1504918807, 783551873, 3082640443, 3233442989, 3988292384, 2596254646, 62317068, 1957810842, 3939845945, 2647816111, 81470997, 1943803523, 3814918930, 2489596804, 225274430, 2053790376, 3826175755, 2466906013, 167816743, 2097651377, 4027552580, 2265490386, 503444072, 1762050814, 4150417245, 2154129355, 426522225, 1852507879, 4275313526, 2312317920, 282753626, 1742555852, 4189708143, 2394877945, 397917763, 1622183637, 3604390888, 2714866558, 953729732, 1340076626, 3518719985, 2797360999, 1068828381, 1219638859, 3624741850, 2936675148, 906185462, 1090812512, 3747672003, 2825379669, 829329135, 1181335161, 3412177804, 3160834842, 628085408, 1382605366, 3423369109, 3138078467, 570562233, 1426400815, 3317316542, 2998733608, 733239954, 1555261956, 3268935591, 3050360625, 752459403, 1541320221, 2607071920, 3965973030, 1969922972, 40735498, 2617837225, 3943577151, 1913087877, 83908371, 2512341634, 3803740692, 2075208622, 213261112, 2463272603, 3855990285, 2094854071, 198958881, 2262029012, 4057260610, 1759359992, 534414190, 2176718541, 4139329115, 1873836001, 414664567, 2282248934, 4279200368, 1711684554, 285281116, 2405801727, 4167216745, 1634467795, 376229701, 2685067896, 3608007406, 1308918612, 956543938, 2808555105, 3495958263, 1231636301, 1047427035, 2932959818, 3654703836, 1088359270, 936918e3, 2847714899, 3736837829, 1202900863, 817233897, 3183342108, 3401237130, 1404277552, 615818150, 3134207493, 3453421203, 1423857449, 601450431, 3009837614, 3294710456, 1567103746, 711928724, 3020668471, 3272380065, 1510334235, 755167117];
        var crc = -1,
            y = 0,
            n = a.length,
            i;
        for (i = from; i < to; i++) {
            y = (crc ^ a[i]) & 255;
            crc = crc >>> 8 ^ table[y]
        }
        return crc ^ -1
    }
    var h = img[0].length,
        w = img[0][0].length,
        s1, s2, next, k, length, a, b, i, j, adler32, crc32;
    var stream = [137, 80, 78, 71, 13, 10, 26, 10, 0, 0, 0, 13, 73, 72, 68, 82, w >> 24 & 255, w >> 16 & 255, w >> 8 & 255, w & 255, h >> 24 & 255, h >> 16 & 255, h >> 8 & 255, h & 255, 8, 2, 0, 0, 0, -1, -2, -3, -4, -5, -6, -7, -8, 73, 68, 65, 84, 8, 29];
    crc32 = crc32Array(stream, 12, 29);
    stream[29] = crc32 >> 24 & 255;
    stream[30] = crc32 >> 16 & 255;
    stream[31] = crc32 >> 8 & 255;
    stream[32] = crc32 & 255;
    s1 = 1;
    s2 = 0;
    for (i = 0; i < h; i++) {
        if (i < h - 1) {
            stream.push(0)
        } else {
            stream.push(1)
        }
        a = 3 * w + 1 + (i === 0) & 255;
        b = 3 * w + 1 + (i === 0) >> 8 & 255;
        stream.push(a);
        stream.push(b);
        stream.push(~a & 255);
        stream.push(~b & 255);
        if (i === 0) stream.push(0);
        for (j = 0; j < w; j++) {
            for (k = 0; k < 3; k++) {
                a = img[k][i][j];
                if (a > 255) a = 255;
                else if (a < 0) a = 0;
                else a = Math.round(a);
                s1 = (s1 + a) % 65521;
                s2 = (s2 + s1) % 65521;
                stream.push(a)
            }
        }
        stream.push(0)
    }
    adler32 = (s2 << 16) + s1;
    stream.push(adler32 >> 24 & 255);
    stream.push(adler32 >> 16 & 255);
    stream.push(adler32 >> 8 & 255);
    stream.push(adler32 & 255);
    length = stream.length - 41;
    stream[33] = length >> 24 & 255;
    stream[34] = length >> 16 & 255;
    stream[35] = length >> 8 & 255;
    stream[36] = length & 255;
    crc32 = crc32Array(stream, 37);
    stream.push(crc32 >> 24 & 255);
    stream.push(crc32 >> 16 & 255);
    stream.push(crc32 >> 8 & 255);
    stream.push(crc32 & 255);
    stream.push(0);
    stream.push(0);
    stream.push(0);
    stream.push(0);
    stream.push(73);
    stream.push(69);
    stream.push(78);
    stream.push(68);
    stream.push(174);
    stream.push(66);
    stream.push(96);
    stream.push(130);
    return "data:image/png;base64," + base64(stream)
};
numeric._dim = function _dim(x) {
    var ret = [];
    while (typeof x === "object") {
        ret.push(x.length);
        x = x[0]
    }
    return ret
};
numeric.dim = function dim(x) {
    var y, z;
    if (typeof x === "object") {
        y = x[0];
        if (typeof y === "object") {
            z = y[0];
            if (typeof z === "object") {
                return numeric._dim(x)
            }
            return [x.length, y.length]
        }
        return [x.length]
    }
    return []
};
numeric.mapreduce = function mapreduce(body, init) {
    return Function("x", "accum", "_s", "_k", 'if(typeof accum === "undefined") accum = ' + init + ";\n" + 'if(typeof x === "number") { var xi = x; ' + body + "; return accum; }\n" + 'if(typeof _s === "undefined") _s = numeric.dim(x);\n' + 'if(typeof _k === "undefined") _k = 0;\n' + "var _n = _s[_k];\n" + "var i,xi;\n" + "if(_k < _s.length-1) {\n" + "    for(i=_n-1;i>=0;i--) {\n" + "        accum = arguments.callee(x[i],accum,_s,_k+1);\n" + "    }" + "    return accum;\n" + "}\n" + "for(i=_n-1;i>=1;i-=2) { \n" + "    xi = x[i];\n" + "    " + body + ";\n" + "    xi = x[i-1];\n" + "    " + body + ";\n" + "}\n" + "if(i === 0) {\n" + "    xi = x[i];\n" + "    " + body + "\n" + "}\n" + "return accum;")
};
numeric.mapreduce2 = function mapreduce2(body, setup) {
    return Function("x", "var n = x.length;\n" + "var i,xi;\n" + setup + ";\n" + "for(i=n-1;i!==-1;--i) { \n" + "    xi = x[i];\n" + "    " + body + ";\n" + "}\n" + "return accum;")
};
numeric.same = function same(x, y) {
    var i, n;
    if (!(x instanceof Array) || !(y instanceof Array)) {
        return false
    }
    n = x.length;
    if (n !== y.length) {
        return false
    }
    for (i = 0; i < n; i++) {
        if (x[i] === y[i]) {
            continue
        }
        if (typeof x[i] === "object") {
            if (!same(x[i], y[i])) return false
        } else {
            return false
        }
    }
    return true
};
numeric.rep = function rep(s, v, k) {
    if (typeof k === "undefined") {
        k = 0
    }
    var n = s[k],
        ret = Array(n),
        i;
    if (k === s.length - 1) {
        for (i = n - 2; i >= 0; i -= 2) {
            ret[i + 1] = v;
            ret[i] = v
        }
        if (i === -1) {
            ret[0] = v
        }
        return ret
    }
    for (i = n - 1; i >= 0; i--) {
        ret[i] = numeric.rep(s, v, k + 1)
    }
    return ret
};
numeric.dotMMsmall = function dotMMsmall(x, y) {
    var i, j, k, p, q, r, ret, foo, bar, woo, i0, k0, p0, r0;
    p = x.length;
    q = y.length;
    r = y[0].length;
    ret = Array(p);
    for (i = p - 1; i >= 0; i--) {
        foo = Array(r);
        bar = x[i];
        for (k = r - 1; k >= 0; k--) {
            woo = bar[q - 1] * y[q - 1][k];
            for (j = q - 2; j >= 1; j -= 2) {
                i0 = j - 1;
                woo += bar[j] * y[j][k] + bar[i0] * y[i0][k]
            }
            if (j === 0) {
                woo += bar[0] * y[0][k]
            }
            foo[k] = woo
        }
        ret[i] = foo
    }
    return ret
};
numeric._getCol = function _getCol(A, j, x) {
    var n = A.length,
        i;
    for (i = n - 1; i > 0; --i) {
        x[i] = A[i][j];
        --i;
        x[i] = A[i][j]
    }
    if (i === 0) x[0] = A[0][j]
};
numeric.dotMMbig = function dotMMbig(x, y) {
    var gc = numeric._getCol,
        p = y.length,
        v = Array(p);
    var m = x.length,
        n = y[0].length,
        A = new Array(m),
        xj;
    var VV = numeric.dotVV;
    var i, j, k, z;
    --p;
    --m;
    for (i = m; i !== -1; --i) A[i] = Array(n);
    --n;
    for (i = n; i !== -1; --i) {
        gc(y, i, v);
        for (j = m; j !== -1; --j) {
            z = 0;
            xj = x[j];
            A[j][i] = VV(xj, v)
        }
    }
    return A
};
numeric.dotMV = function dotMV(x, y) {
    var p = x.length,
        q = y.length,
        i;
    var ret = Array(p),
        dotVV = numeric.dotVV;
    for (i = p - 1; i >= 0; i--) {
        ret[i] = dotVV(x[i], y)
    }
    return ret
};
numeric.dotVM = function dotVM(x, y) {
    var i, j, k, p, q, r, ret, foo, bar, woo, i0, k0, p0, r0, s1, s2, s3, baz, accum;
    p = x.length;
    q = y[0].length;
    ret = Array(q);
    for (k = q - 1; k >= 0; k--) {
        woo = x[p - 1] * y[p - 1][k];
        for (j = p - 2; j >= 1; j -= 2) {
            i0 = j - 1;
            woo += x[j] * y[j][k] + x[i0] * y[i0][k]
        }
        if (j === 0) {
            woo += x[0] * y[0][k]
        }
        ret[k] = woo
    }
    return ret
};
numeric.dotVV = function dotVV(x, y) {
    var i, n = x.length,
        i1, ret = x[n - 1] * y[n - 1];
    for (i = n - 2; i >= 1; i -= 2) {
        i1 = i - 1;
        ret += x[i] * y[i] + x[i1] * y[i1]
    }
    if (i === 0) {
        ret += x[0] * y[0]
    }
    return ret
};
numeric.dot = function dot(x, y) {
    var d = numeric.dim;
    switch (d(x).length * 1e3 + d(y).length) {
        case 2002:
            if (y.length < 10) return numeric.dotMMsmall(x, y);
            else return numeric.dotMMbig(x, y);
        case 2001:
            return numeric.dotMV(x, y);
        case 1002:
            return numeric.dotVM(x, y);
        case 1001:
            return numeric.dotVV(x, y);
        case 1e3:
            return numeric.mulVS(x, y);
        case 1:
            return numeric.mulSV(x, y);
        case 0:
            return x * y;
        default:
            throw new Error("numeric.dot only works on vectors and matrices")
    }
};
numeric.diag = function diag(d) {
    var i, i1, j, n = d.length,
        A = Array(n),
        Ai;
    for (i = n - 1; i >= 0; i--) {
        Ai = Array(n);
        i1 = i + 2;
        for (j = n - 1; j >= i1; j -= 2) {
            Ai[j] = 0;
            Ai[j - 1] = 0
        }
        if (j > i) {
            Ai[j] = 0
        }
        Ai[i] = d[i];
        for (j = i - 1; j >= 1; j -= 2) {
            Ai[j] = 0;
            Ai[j - 1] = 0
        }
        if (j === 0) {
            Ai[0] = 0
        }
        A[i] = Ai
    }
    return A
};
numeric.getDiag = function(A) {
    var n = Math.min(A.length, A[0].length),
        i, ret = Array(n);
    for (i = n - 1; i >= 1; --i) {
        ret[i] = A[i][i];
        --i;
        ret[i] = A[i][i]
    }
    if (i === 0) {
        ret[0] = A[0][0]
    }
    return ret
};
numeric.identity = function identity(n) {
    return numeric.diag(numeric.rep([n], 1))
};
numeric.pointwise = function pointwise(params, body, setup) {
    if (typeof setup === "undefined") {
        setup = ""
    }
    var fun = [];
    var k;
    var avec = /\[i\]$/,
        p, thevec = "";
    var haveret = false;
    for (k = 0; k < params.length; k++) {
        if (avec.test(params[k])) {
            p = params[k].substring(0, params[k].length - 3);
            thevec = p
        } else {
            p = params[k]
        }
        if (p === "ret") haveret = true;
        fun.push(p)
    }
    fun[params.length] = "_s";
    fun[params.length + 1] = "_k";
    fun[params.length + 2] = 'if(typeof _s === "undefined") _s = numeric.dim(' + thevec + ");\n" + 'if(typeof _k === "undefined") _k = 0;\n' + "var _n = _s[_k];\n" + "var i" + (haveret ? "" : ", ret = Array(_n)") + ";\n" + "if(_k < _s.length-1) {\n" + "    for(i=_n-1;i>=0;i--) ret[i] = arguments.callee(" + params.join(",") + ",_s,_k+1);\n" + "    return ret;\n" + "}\n" + setup + "\n" + "for(i=_n-1;i!==-1;--i) {\n" + "    " + body + "\n" + "}\n" + "return ret;";
    return Function.apply(null, fun)
};
numeric.pointwise2 = function pointwise2(params, body, setup) {
    if (typeof setup === "undefined") {
        setup = ""
    }
    var fun = [];
    var k;
    var avec = /\[i\]$/,
        p, thevec = "";
    var haveret = false;
    for (k = 0; k < params.length; k++) {
        if (avec.test(params[k])) {
            p = params[k].substring(0, params[k].length - 3);
            thevec = p
        } else {
            p = params[k]
        }
        if (p === "ret") haveret = true;
        fun.push(p)
    }
    fun[params.length] = "var _n = " + thevec + ".length;\n" + "var i" + (haveret ? "" : ", ret = Array(_n)") + ";\n" + setup + "\n" + "for(i=_n-1;i!==-1;--i) {\n" + body + "\n" + "}\n" + "return ret;";
    return Function.apply(null, fun)
};
numeric._biforeach = function _biforeach(x, y, s, k, f) {
    if (k === s.length - 1) {
        f(x, y);
        return
    }
    var i, n = s[k];
    for (i = n - 1; i >= 0; i--) {
        _biforeach(typeof x === "object" ? x[i] : x, typeof y === "object" ? y[i] : y, s, k + 1, f)
    }
};
numeric._biforeach2 = function _biforeach2(x, y, s, k, f) {
    if (k === s.length - 1) {
        return f(x, y)
    }
    var i, n = s[k],
        ret = Array(n);
    for (i = n - 1; i >= 0; --i) {
        ret[i] = _biforeach2(typeof x === "object" ? x[i] : x, typeof y === "object" ? y[i] : y, s, k + 1, f)
    }
    return ret
};
numeric._foreach = function _foreach(x, s, k, f) {
    if (k === s.length - 1) {
        f(x);
        return
    }
    var i, n = s[k];
    for (i = n - 1; i >= 0; i--) {
        _foreach(x[i], s, k + 1, f)
    }
};
numeric._foreach2 = function _foreach2(x, s, k, f) {
    if (k === s.length - 1) {
        return f(x)
    }
    var i, n = s[k],
        ret = Array(n);
    for (i = n - 1; i >= 0; i--) {
        ret[i] = _foreach2(x[i], s, k + 1, f)
    }
    return ret
};
numeric.ops2 = {
    add: "+",
    sub: "-",
    mul: "*",
    div: "/",
    mod: "%",
    and: "&&",
    or: "||",
    eq: "===",
    neq: "!==",
    lt: "<",
    gt: ">",
    leq: "<=",
    geq: ">=",
    band: "&",
    bor: "|",
    bxor: "^",
    lshift: "<<",
    rshift: ">>",
    rrshift: ">>>"
};
numeric.opseq = {
    addeq: "+=",
    subeq: "-=",
    muleq: "*=",
    diveq: "/=",
    modeq: "%=",
    lshifteq: "<<=",
    rshifteq: ">>=",
    rrshifteq: ">>>=",
    bandeq: "&=",
    boreq: "|=",
    bxoreq: "^="
};
numeric.mathfuns = ["abs", "acos", "asin", "atan", "ceil", "cos", "exp", "floor", "log", "round", "sin", "sqrt", "tan", "isNaN", "isFinite"];
numeric.mathfuns2 = ["atan2", "pow", "max", "min"];
numeric.ops1 = {
    neg: "-",
    not: "!",
    bnot: "~",
    clone: ""
};
numeric.mapreducers = {
    any: ["if(xi) return true;", "var accum = false;"],
    all: ["if(!xi) return false;", "var accum = true;"],
    sum: ["accum += xi;", "var accum = 0;"],
    prod: ["accum *= xi;", "var accum = 1;"],
    norm2Squared: ["accum += xi*xi;", "var accum = 0;"],
    norminf: ["accum = max(accum,abs(xi));", "var accum = 0, max = Math.max, abs = Math.abs;"],
    norm1: ["accum += abs(xi)", "var accum = 0, abs = Math.abs;"],
    sup: ["accum = max(accum,xi);", "var accum = -Infinity, max = Math.max;"],
    inf: ["accum = min(accum,xi);", "var accum = Infinity, min = Math.min;"]
};
(function() {
    var i, o;
    for (i = 0; i < numeric.mathfuns2.length; ++i) {
        o = numeric.mathfuns2[i];
        numeric.ops2[o] = o
    }
    for (i in numeric.ops2) {
        if (numeric.ops2.hasOwnProperty(i)) {
            o = numeric.ops2[i];
            var code, codeeq, setup = "";
            if (numeric.myIndexOf.call(numeric.mathfuns2, i) !== -1) {
                setup = "var " + o + " = Math." + o + ";\n";
                code = function(r, x, y) {
                    return r + " = " + o + "(" + x + "," + y + ")"
                };
                codeeq = function(x, y) {
                    return x + " = " + o + "(" + x + "," + y + ")"
                }
            } else {
                code = function(r, x, y) {
                    return r + " = " + x + " " + o + " " + y
                };
                if (numeric.opseq.hasOwnProperty(i + "eq")) {
                    codeeq = function(x, y) {
                        return x + " " + o + "= " + y
                    }
                } else {
                    codeeq = function(x, y) {
                        return x + " = " + x + " " + o + " " + y
                    }
                }
            }
            numeric[i + "VV"] = numeric.pointwise2(["x[i]", "y[i]"], code("ret[i]", "x[i]", "y[i]"), setup);
            numeric[i + "SV"] = numeric.pointwise2(["x", "y[i]"], code("ret[i]", "x", "y[i]"), setup);
            numeric[i + "VS"] = numeric.pointwise2(["x[i]", "y"], code("ret[i]", "x[i]", "y"), setup);
            numeric[i] = Function("var n = arguments.length, i, x = arguments[0], y;\n" + "var VV = numeric." + i + "VV, VS = numeric." + i + "VS, SV = numeric." + i + "SV;\n" + "var dim = numeric.dim;\n" + "for(i=1;i!==n;++i) { \n" + "  y = arguments[i];\n" + '  if(typeof x === "object") {\n' + '      if(typeof y === "object") x = numeric._biforeach2(x,y,dim(x),0,VV);\n' + "      else x = numeric._biforeach2(x,y,dim(x),0,VS);\n" + '  } else if(typeof y === "object") x = numeric._biforeach2(x,y,dim(y),0,SV);\n' + "  else " + codeeq("x", "y") + "\n" + "}\nreturn x;\n");
            numeric[o] = numeric[i];
            numeric[i + "eqV"] = numeric.pointwise2(["ret[i]", "x[i]"], codeeq("ret[i]", "x[i]"), setup);
            numeric[i + "eqS"] = numeric.pointwise2(["ret[i]", "x"], codeeq("ret[i]", "x"), setup);
            numeric[i + "eq"] = Function("var n = arguments.length, i, x = arguments[0], y;\n" + "var V = numeric." + i + "eqV, S = numeric." + i + "eqS\n" + "var s = numeric.dim(x);\n" + "for(i=1;i!==n;++i) { \n" + "  y = arguments[i];\n" + '  if(typeof y === "object") numeric._biforeach(x,y,s,0,V);\n' + "  else numeric._biforeach(x,y,s,0,S);\n" + "}\nreturn x;\n")
        }
    }
    for (i = 0; i < numeric.mathfuns2.length; ++i) {
        o = numeric.mathfuns2[i];
        delete numeric.ops2[o]
    }
    for (i = 0; i < numeric.mathfuns.length; ++i) {
        o = numeric.mathfuns[i];
        numeric.ops1[o] = o
    }
    for (i in numeric.ops1) {
        if (numeric.ops1.hasOwnProperty(i)) {
            setup = "";
            o = numeric.ops1[i];
            if (numeric.myIndexOf.call(numeric.mathfuns, i) !== -1) {
                if (Math.hasOwnProperty(o)) setup = "var " + o + " = Math." + o + ";\n"
            }
            numeric[i + "eqV"] = numeric.pointwise2(["ret[i]"], "ret[i] = " + o + "(ret[i]);", setup);
            numeric[i + "eq"] = Function("x", 'if(typeof x !== "object") return ' + o + "x\n" + "var i;\n" + "var V = numeric." + i + "eqV;\n" + "var s = numeric.dim(x);\n" + "numeric._foreach(x,s,0,V);\n" + "return x;\n");
            numeric[i + "V"] = numeric.pointwise2(["x[i]"], "ret[i] = " + o + "(x[i]);", setup);
            numeric[i] = Function("x", 'if(typeof x !== "object") return ' + o + "(x)\n" + "var i;\n" + "var V = numeric." + i + "V;\n" + "var s = numeric.dim(x);\n" + "return numeric._foreach2(x,s,0,V);\n")
        }
    }
    for (i = 0; i < numeric.mathfuns.length; ++i) {
        o = numeric.mathfuns[i];
        delete numeric.ops1[o]
    }
    for (i in numeric.mapreducers) {
        if (numeric.mapreducers.hasOwnProperty(i)) {
            o = numeric.mapreducers[i];
            numeric[i + "V"] = numeric.mapreduce2(o[0], o[1]);
            numeric[i] = Function("x", "s", "k", o[1] + 'if(typeof x !== "object") {' + "    xi = x;\n" + o[0] + ";\n" + "    return accum;\n" + "}" + 'if(typeof s === "undefined") s = numeric.dim(x);\n' + 'if(typeof k === "undefined") k = 0;\n' + "if(k === s.length-1) return numeric." + i + "V(x);\n" + "var xi;\n" + "var n = x.length, i;\n" + "for(i=n-1;i!==-1;--i) {\n" + "   xi = arguments.callee(x[i]);\n" + o[0] + ";\n" + "}\n" + "return accum;\n")
        }
    }
})();
numeric.truncVV = numeric.pointwise(["x[i]", "y[i]"], "ret[i] = round(x[i]/y[i])*y[i];", "var round = Math.round;");
numeric.truncVS = numeric.pointwise(["x[i]", "y"], "ret[i] = round(x[i]/y)*y;", "var round = Math.round;");
numeric.truncSV = numeric.pointwise(["x", "y[i]"], "ret[i] = round(x/y[i])*y[i];", "var round = Math.round;");
numeric.trunc = function trunc(x, y) {
    if (typeof x === "object") {
        if (typeof y === "object") return numeric.truncVV(x, y);
        return numeric.truncVS(x, y)
    }
    if (typeof y === "object") return numeric.truncSV(x, y);
    return Math.round(x / y) * y
};
numeric.inv = function inv(x) {
    var s = numeric.dim(x),
        abs = Math.abs,
        m = s[0],
        n = s[1];
    var A = numeric.clone(x),
        Ai, Aj;
    var I = numeric.identity(m),
        Ii, Ij;
    var i, j, k, x;
    for (j = 0; j < n; ++j) {
        var i0 = -1;
        var v0 = -1;
        for (i = j; i !== m; ++i) {
            k = abs(A[i][j]);
            if (k > v0) {
                i0 = i;
                v0 = k
            }
        }
        Aj = A[i0];
        A[i0] = A[j];
        A[j] = Aj;
        Ij = I[i0];
        I[i0] = I[j];
        I[j] = Ij;
        x = Aj[j];
        for (k = j; k !== n; ++k) Aj[k] /= x;
        for (k = n - 1; k !== -1; --k) Ij[k] /= x;
        for (i = m - 1; i !== -1; --i) {
            if (i !== j) {
                Ai = A[i];
                Ii = I[i];
                x = Ai[j];
                for (k = j + 1; k !== n; ++k) Ai[k] -= Aj[k] * x;
                for (k = n - 1; k > 0; --k) {
                    Ii[k] -= Ij[k] * x;
                    --k;
                    Ii[k] -= Ij[k] * x
                }
                if (k === 0) Ii[0] -= Ij[0] * x
            }
        }
    }
    return I
};
numeric.det = function det(x) {
    var s = numeric.dim(x);
    if (s.length !== 2 || s[0] !== s[1]) {
        throw new Error("numeric: det() only works on square matrices")
    }
    var n = s[0],
        ret = 1,
        i, j, k, A = numeric.clone(x),
        Aj, Ai, alpha, temp, k1, k2, k3;
    for (j = 0; j < n - 1; j++) {
        k = j;
        for (i = j + 1; i < n; i++) {
            if (Math.abs(A[i][j]) > Math.abs(A[k][j])) {
                k = i
            }
        }
        if (k !== j) {
            temp = A[k];
            A[k] = A[j];
            A[j] = temp;
            ret *= -1
        }
        Aj = A[j];
        for (i = j + 1; i < n; i++) {
            Ai = A[i];
            alpha = Ai[j] / Aj[j];
            for (k = j + 1; k < n - 1; k += 2) {
                k1 = k + 1;
                Ai[k] -= Aj[k] * alpha;
                Ai[k1] -= Aj[k1] * alpha
            }
            if (k !== n) {
                Ai[k] -= Aj[k] * alpha
            }
        }
        if (Aj[j] === 0) {
            return 0
        }
        ret *= Aj[j]
    }
    return ret * A[j][j]
};
numeric.transpose = function transpose(x) {
    var i, j, m = x.length,
        n = x[0].length,
        ret = Array(n),
        A0, A1, Bj;
    for (j = 0; j < n; j++) ret[j] = Array(m);
    for (i = m - 1; i >= 1; i -= 2) {
        A1 = x[i];
        A0 = x[i - 1];
        for (j = n - 1; j >= 1; --j) {
            Bj = ret[j];
            Bj[i] = A1[j];
            Bj[i - 1] = A0[j];
            --j;
            Bj = ret[j];
            Bj[i] = A1[j];
            Bj[i - 1] = A0[j]
        }
        if (j === 0) {
            Bj = ret[0];
            Bj[i] = A1[0];
            Bj[i - 1] = A0[0]
        }
    }
    if (i === 0) {
        A0 = x[0];
        for (j = n - 1; j >= 1; --j) {
            ret[j][0] = A0[j];
            --j;
            ret[j][0] = A0[j]
        }
        if (j === 0) {
            ret[0][0] = A0[0]
        }
    }
    return ret
};
numeric.negtranspose = function negtranspose(x) {
    var i, j, m = x.length,
        n = x[0].length,
        ret = Array(n),
        A0, A1, Bj;
    for (j = 0; j < n; j++) ret[j] = Array(m);
    for (i = m - 1; i >= 1; i -= 2) {
        A1 = x[i];
        A0 = x[i - 1];
        for (j = n - 1; j >= 1; --j) {
            Bj = ret[j];
            Bj[i] = -A1[j];
            Bj[i - 1] = -A0[j];
            --j;
            Bj = ret[j];
            Bj[i] = -A1[j];
            Bj[i - 1] = -A0[j]
        }
        if (j === 0) {
            Bj = ret[0];
            Bj[i] = -A1[0];
            Bj[i - 1] = -A0[0]
        }
    }
    if (i === 0) {
        A0 = x[0];
        for (j = n - 1; j >= 1; --j) {
            ret[j][0] = -A0[j];
            --j;
            ret[j][0] = -A0[j]
        }
        if (j === 0) {
            ret[0][0] = -A0[0]
        }
    }
    return ret
};
numeric._random = function _random(s, k) {
    var i, n = s[k],
        ret = Array(n),
        rnd;
    if (k === s.length - 1) {
        rnd = Math.random;
        for (i = n - 1; i >= 1; i -= 2) {
            ret[i] = rnd();
            ret[i - 1] = rnd()
        }
        if (i === 0) {
            ret[0] = rnd()
        }
        return ret
    }
    for (i = n - 1; i >= 0; i--) ret[i] = _random(s, k + 1);
    return ret
};
numeric.random = function random(s) {
    return numeric._random(s, 0)
};
numeric.norm2 = function norm2(x) {
    return Math.sqrt(numeric.norm2Squared(x))
};
numeric.linspace = function linspace(a, b, n) {
    if (typeof n === "undefined") n = Math.max(Math.round(b - a) + 1, 1);
    if (n < 2) {
        return n === 1 ? [a] : []
    }
    var i, ret = Array(n);
    n--;
    for (i = n; i >= 0; i--) {
        ret[i] = (i * b + (n - i) * a) / n
    }
    return ret
};
numeric.getBlock = function getBlock(x, from, to) {
    var s = numeric.dim(x);

    function foo(x, k) {
        var i, a = from[k],
            n = to[k] - a,
            ret = Array(n);
        if (k === s.length - 1) {
            for (i = n; i >= 0; i--) {
                ret[i] = x[i + a]
            }
            return ret
        }
        for (i = n; i >= 0; i--) {
            ret[i] = foo(x[i + a], k + 1)
        }
        return ret
    }
    return foo(x, 0)
};
numeric.setBlock = function setBlock(x, from, to, B) {
    var s = numeric.dim(x);

    function foo(x, y, k) {
        var i, a = from[k],
            n = to[k] - a;
        if (k === s.length - 1) {
            for (i = n; i >= 0; i--) {
                x[i + a] = y[i]
            }
        }
        for (i = n; i >= 0; i--) {
            foo(x[i + a], y[i], k + 1)
        }
    }
    foo(x, B, 0);
    return x
};
numeric.getRange = function getRange(A, I, J) {
    var m = I.length,
        n = J.length;
    var i, j;
    var B = Array(m),
        Bi, AI;
    for (i = m - 1; i !== -1; --i) {
        B[i] = Array(n);
        Bi = B[i];
        AI = A[I[i]];
        for (j = n - 1; j !== -1; --j) Bi[j] = AI[J[j]]
    }
    return B
};
numeric.blockMatrix = function blockMatrix(X) {
    var s = numeric.dim(X);
    if (s.length < 4) return numeric.blockMatrix([X]);
    var m = s[0],
        n = s[1],
        M, N, i, j, Xij;
    M = 0;
    N = 0;
    for (i = 0; i < m; ++i) M += X[i][0].length;
    for (j = 0; j < n; ++j) N += X[0][j][0].length;
    var Z = Array(M);
    for (i = 0; i < M; ++i) Z[i] = Array(N);
    var I = 0,
        J, ZI, k, l, Xijk;
    for (i = 0; i < m; ++i) {
        J = N;
        for (j = n - 1; j !== -1; --j) {
            Xij = X[i][j];
            J -= Xij[0].length;
            for (k = Xij.length - 1; k !== -1; --k) {
                Xijk = Xij[k];
                ZI = Z[I + k];
                for (l = Xijk.length - 1; l !== -1; --l) ZI[J + l] = Xijk[l]
            }
        }
        I += X[i][0].length
    }
    return Z
};
numeric.tensor = function tensor(x, y) {
    if (typeof x === "number" || typeof y === "number") return numeric.mul(x, y);
    var s1 = numeric.dim(x),
        s2 = numeric.dim(y);
    if (s1.length !== 1 || s2.length !== 1) {
        throw new Error("numeric: tensor product is only defined for vectors")
    }
    var m = s1[0],
        n = s2[0],
        A = Array(m),
        Ai, i, j, xi;
    for (i = m - 1; i >= 0; i--) {
        Ai = Array(n);
        xi = x[i];
        for (j = n - 1; j >= 3; --j) {
            Ai[j] = xi * y[j];
            --j;
            Ai[j] = xi * y[j];
            --j;
            Ai[j] = xi * y[j];
            --j;
            Ai[j] = xi * y[j]
        }
        while (j >= 0) {
            Ai[j] = xi * y[j];
            --j
        }
        A[i] = Ai
    }
    return A
};
numeric.T = function T(x, y) {
    this.x = x;
    this.y = y
};
numeric.t = function t(x, y) {
    return new numeric.T(x, y)
};
numeric.Tbinop = function Tbinop(rr, rc, cr, cc, setup) {
    var io = numeric.indexOf;
    if (typeof setup !== "string") {
        var k;
        setup = "";
        for (k in numeric) {
            if (numeric.hasOwnProperty(k) && (rr.indexOf(k) >= 0 || rc.indexOf(k) >= 0 || cr.indexOf(k) >= 0 || cc.indexOf(k) >= 0) && k.length > 1) {
                setup += "var " + k + " = numeric." + k + ";\n"
            }
        }
    }
    return Function(["y"], "var x = this;\n" + "if(!(y instanceof numeric.T)) { y = new numeric.T(y); }\n" + setup + "\n" + "if(x.y) {" + "  if(y.y) {" + "    return new numeric.T(" + cc + ");\n" + "  }\n" + "  return new numeric.T(" + cr + ");\n" + "}\n" + "if(y.y) {\n" + "  return new numeric.T(" + rc + ");\n" + "}\n" + "return new numeric.T(" + rr + ");\n")
};
numeric.T.prototype.add = numeric.Tbinop("add(x.x,y.x)", "add(x.x,y.x),y.y", "add(x.x,y.x),x.y", "add(x.x,y.x),add(x.y,y.y)");
numeric.T.prototype.sub = numeric.Tbinop("sub(x.x,y.x)", "sub(x.x,y.x),neg(y.y)", "sub(x.x,y.x),x.y", "sub(x.x,y.x),sub(x.y,y.y)");
numeric.T.prototype.mul = numeric.Tbinop("mul(x.x,y.x)", "mul(x.x,y.x),mul(x.x,y.y)", "mul(x.x,y.x),mul(x.y,y.x)", "sub(mul(x.x,y.x),mul(x.y,y.y)),add(mul(x.x,y.y),mul(x.y,y.x))");
numeric.T.prototype.reciprocal = function reciprocal() {
    var mul = numeric.mul,
        div = numeric.div;
    if (this.y) {
        var d = numeric.add(mul(this.x, this.x), mul(this.y, this.y));
        return new numeric.T(div(this.x, d), div(numeric.neg(this.y), d))
    }
    return new T(div(1, this.x))
};
numeric.T.prototype.div = function div(y) {
    if (!(y instanceof numeric.T)) y = new numeric.T(y);
    if (y.y) {
        return this.mul(y.reciprocal())
    }
    var div = numeric.div;
    if (this.y) {
        return new numeric.T(div(this.x, y.x), div(this.y, y.x))
    }
    return new numeric.T(div(this.x, y.x))
};
numeric.T.prototype.dot = numeric.Tbinop("dot(x.x,y.x)", "dot(x.x,y.x),dot(x.x,y.y)", "dot(x.x,y.x),dot(x.y,y.x)", "sub(dot(x.x,y.x),dot(x.y,y.y)),add(dot(x.x,y.y),dot(x.y,y.x))");
numeric.T.prototype.transpose = function transpose() {
    var t = numeric.transpose,
        x = this.x,
        y = this.y;
    if (y) {
        return new numeric.T(t(x), t(y))
    }
    return new numeric.T(t(x))
};
numeric.T.prototype.transjugate = function transjugate() {
    var t = numeric.transpose,
        x = this.x,
        y = this.y;
    if (y) {
        return new numeric.T(t(x), numeric.negtranspose(y))
    }
    return new numeric.T(t(x))
};
numeric.Tunop = function Tunop(r, c, s) {
    if (typeof s !== "string") {
        s = ""
    }
    return Function("var x = this;\n" + s + "\n" + "if(x.y) {" + "  " + c + ";\n" + "}\n" + r + ";\n")
};
numeric.T.prototype.exp = numeric.Tunop("return new numeric.T(ex)", "return new numeric.T(mul(cos(x.y),ex),mul(sin(x.y),ex))", "var ex = numeric.exp(x.x), cos = numeric.cos, sin = numeric.sin, mul = numeric.mul;");
numeric.T.prototype.conj = numeric.Tunop("return new numeric.T(x.x);", "return new numeric.T(x.x,numeric.neg(x.y));");
numeric.T.prototype.neg = numeric.Tunop("return new numeric.T(neg(x.x));", "return new numeric.T(neg(x.x),neg(x.y));", "var neg = numeric.neg;");
numeric.T.prototype.sin = numeric.Tunop("return new numeric.T(numeric.sin(x.x))", "return x.exp().sub(x.neg().exp()).div(new numeric.T(0,2));");
numeric.T.prototype.cos = numeric.Tunop("return new numeric.T(numeric.cos(x.x))", "return x.exp().add(x.neg().exp()).div(2);");
numeric.T.prototype.abs = numeric.Tunop("return new numeric.T(numeric.abs(x.x));", "return new numeric.T(numeric.sqrt(numeric.add(mul(x.x,x.x),mul(x.y,x.y))));", "var mul = numeric.mul;");
numeric.T.prototype.log = numeric.Tunop("return new numeric.T(numeric.log(x.x));", "var theta = new numeric.T(numeric.atan2(x.y,x.x)), r = x.abs();\n" + "return new numeric.T(numeric.log(r.x),theta.x);");
numeric.T.prototype.norm2 = numeric.Tunop("return numeric.norm2(x.x);", "var f = numeric.norm2Squared;\n" + "return Math.sqrt(f(x.x)+f(x.y));");
numeric.T.prototype.inv = function inv() {
    var A = this;
    if (typeof A.y === "undefined") {
        return new numeric.T(numeric.inv(A.x))
    }
    var n = A.x.length,
        i, j, k;
    var Rx = numeric.identity(n),
        Ry = numeric.rep([n, n], 0);
    var Ax = numeric.clone(A.x),
        Ay = numeric.clone(A.y);
    var Aix, Aiy, Ajx, Ajy, Rix, Riy, Rjx, Rjy;
    var i, j, k, d, d1, ax, ay, bx, by, temp;
    for (i = 0; i < n; i++) {
        ax = Ax[i][i];
        ay = Ay[i][i];
        d = ax * ax + ay * ay;
        k = i;
        for (j = i + 1; j < n; j++) {
            ax = Ax[j][i];
            ay = Ay[j][i];
            d1 = ax * ax + ay * ay;
            if (d1 > d) {
                k = j;
                d = d1
            }
        }
        if (k !== i) {
            temp = Ax[i];
            Ax[i] = Ax[k];
            Ax[k] = temp;
            temp = Ay[i];
            Ay[i] = Ay[k];
            Ay[k] = temp;
            temp = Rx[i];
            Rx[i] = Rx[k];
            Rx[k] = temp;
            temp = Ry[i];
            Ry[i] = Ry[k];
            Ry[k] = temp
        }
        Aix = Ax[i];
        Aiy = Ay[i];
        Rix = Rx[i];
        Riy = Ry[i];
        ax = Aix[i];
        ay = Aiy[i];
        for (j = i + 1; j < n; j++) {
            bx = Aix[j];
            by = Aiy[j];
            Aix[j] = (bx * ax + by * ay) / d;
            Aiy[j] = (by * ax - bx * ay) / d
        }
        for (j = 0; j < n; j++) {
            bx = Rix[j];
            by = Riy[j];
            Rix[j] = (bx * ax + by * ay) / d;
            Riy[j] = (by * ax - bx * ay) / d
        }
        for (j = i + 1; j < n; j++) {
            Ajx = Ax[j];
            Ajy = Ay[j];
            Rjx = Rx[j];
            Rjy = Ry[j];
            ax = Ajx[i];
            ay = Ajy[i];
            for (k = i + 1; k < n; k++) {
                bx = Aix[k];
                by = Aiy[k];
                Ajx[k] -= bx * ax - by * ay;
                Ajy[k] -= by * ax + bx * ay
            }
            for (k = 0; k < n; k++) {
                bx = Rix[k];
                by = Riy[k];
                Rjx[k] -= bx * ax - by * ay;
                Rjy[k] -= by * ax + bx * ay
            }
        }
    }
    for (i = n - 1; i > 0; i--) {
        Rix = Rx[i];
        Riy = Ry[i];
        for (j = i - 1; j >= 0; j--) {
            Rjx = Rx[j];
            Rjy = Ry[j];
            ax = Ax[j][i];
            ay = Ay[j][i];
            for (k = n - 1; k >= 0; k--) {
                bx = Rix[k];
                by = Riy[k];
                Rjx[k] -= ax * bx - ay * by;
                Rjy[k] -= ax * by + ay * bx
            }
        }
    }
    return new numeric.T(Rx, Ry)
};
numeric.T.prototype.get = function get(i) {
    var x = this.x,
        y = this.y,
        k = 0,
        ik, n = i.length;
    if (y) {
        while (k < n) {
            ik = i[k];
            x = x[ik];
            y = y[ik];
            k++
        }
        return new numeric.T(x, y)
    }
    while (k < n) {
        ik = i[k];
        x = x[ik];
        k++
    }
    return new numeric.T(x)
};
numeric.T.prototype.set = function set(i, v) {
    var x = this.x,
        y = this.y,
        k = 0,
        ik, n = i.length,
        vx = v.x,
        vy = v.y;
    if (n === 0) {
        if (vy) {
            this.y = vy
        } else if (y) {
            this.y = undefined
        }
        this.x = x;
        return this
    }
    if (vy) {
        if (y) {} else {
            y = numeric.rep(numeric.dim(x), 0);
            this.y = y
        }
        while (k < n - 1) {
            ik = i[k];
            x = x[ik];
            y = y[ik];
            k++
        }
        ik = i[k];
        x[ik] = vx;
        y[ik] = vy;
        return this
    }
    if (y) {
        while (k < n - 1) {
            ik = i[k];
            x = x[ik];
            y = y[ik];
            k++
        }
        ik = i[k];
        x[ik] = vx;
        if (vx instanceof Array) y[ik] = numeric.rep(numeric.dim(vx), 0);
        else y[ik] = 0;
        return this
    }
    while (k < n - 1) {
        ik = i[k];
        x = x[ik];
        k++
    }
    ik = i[k];
    x[ik] = vx;
    return this
};
numeric.T.prototype.getRows = function getRows(i0, i1) {
    var n = i1 - i0 + 1,
        j;
    var rx = Array(n),
        ry, x = this.x,
        y = this.y;
    for (j = i0; j <= i1; j++) {
        rx[j - i0] = x[j]
    }
    if (y) {
        ry = Array(n);
        for (j = i0; j <= i1; j++) {
            ry[j - i0] = y[j]
        }
        return new numeric.T(rx, ry)
    }
    return new numeric.T(rx)
};
numeric.T.prototype.setRows = function setRows(i0, i1, A) {
    var j;
    var rx = this.x,
        ry = this.y,
        x = A.x,
        y = A.y;
    for (j = i0; j <= i1; j++) {
        rx[j] = x[j - i0]
    }
    if (y) {
        if (!ry) {
            ry = numeric.rep(numeric.dim(rx), 0);
            this.y = ry
        }
        for (j = i0; j <= i1; j++) {
            ry[j] = y[j - i0]
        }
    } else if (ry) {
        for (j = i0; j <= i1; j++) {
            ry[j] = numeric.rep([x[j - i0].length], 0)
        }
    }
    return this
};
numeric.T.prototype.getRow = function getRow(k) {
    var x = this.x,
        y = this.y;
    if (y) {
        return new numeric.T(x[k], y[k])
    }
    return new numeric.T(x[k])
};
numeric.T.prototype.setRow = function setRow(i, v) {
    var rx = this.x,
        ry = this.y,
        x = v.x,
        y = v.y;
    rx[i] = x;
    if (y) {
        if (!ry) {
            ry = numeric.rep(numeric.dim(rx), 0);
            this.y = ry
        }
        ry[i] = y
    } else if (ry) {
        ry = numeric.rep([x.length], 0)
    }
    return this
};
numeric.T.prototype.getBlock = function getBlock(from, to) {
    var x = this.x,
        y = this.y,
        b = numeric.getBlock;
    if (y) {
        return new numeric.T(b(x, from, to), b(y, from, to))
    }
    return new numeric.T(b(x, from, to))
};
numeric.T.prototype.setBlock = function setBlock(from, to, A) {
    if (!(A instanceof numeric.T)) A = new numeric.T(A);
    var x = this.x,
        y = this.y,
        b = numeric.setBlock,
        Ax = A.x,
        Ay = A.y;
    if (Ay) {
        if (!y) {
            this.y = numeric.rep(numeric.dim(this), 0);
            y = this.y
        }
        b(x, from, to, Ax);
        b(y, from, to, Ay);
        return this
    }
    b(x, from, to, Ax);
    if (y) b(y, from, to, numeric.rep(numeric.dim(Ax), 0))
};
numeric.T.rep = function rep(s, v) {
    var T = numeric.T;
    if (!(v instanceof T)) v = new T(v);
    var x = v.x,
        y = v.y,
        r = numeric.rep;
    if (y) return new T(r(s, x), r(s, y));
    return new T(r(s, x))
};
numeric.T.diag = function diag(d) {
    if (!(d instanceof numeric.T)) d = new numeric.T(d);
    var x = d.x,
        y = d.y,
        diag = numeric.diag;
    if (y) return new numeric.T(diag(x), diag(y));
    return new numeric.T(diag(x))
};
numeric.T.eig = function eig() {
    if (this.y) {
        throw new Error("eig: not implemented for complex matrices.")
    }
    return numeric.eig(this.x)
};
numeric.T.identity = function identity(n) {
    return new numeric.T(numeric.identity(n))
};
numeric.T.prototype.getDiag = function getDiag() {
    var n = numeric;
    var x = this.x,
        y = this.y;
    if (y) {
        return new n.T(n.getDiag(x), n.getDiag(y))
    }
    return new n.T(n.getDiag(x))
};
numeric.house = function house(x) {
    var v = numeric.clone(x);
    var s = x[0] >= 0 ? 1 : -1;
    var alpha = s * numeric.norm2(x);
    v[0] += alpha;
    var foo = numeric.norm2(v);
    if (foo === 0) {
        throw new Error("eig: internal error")
    }
    return numeric.div(v, foo)
};
numeric.toUpperHessenberg = function toUpperHessenberg(me) {
    var s = numeric.dim(me);
    if (s.length !== 2 || s[0] !== s[1]) {
        throw new Error("numeric: toUpperHessenberg() only works on square matrices");
    }
    var m = s[0],
        i, j, k, x, v, A = numeric.clone(me),
        B, C, Ai, Ci, Q = numeric.identity(m),
        Qi;
    for (j = 0; j < m - 2; j++) {
        x = Array(m - j - 1);
        for (i = j + 1; i < m; i++) {
            x[i - j - 1] = A[i][j]
        }
        if (numeric.norm2(x) > 0) {
            v = numeric.house(x);
            B = numeric.getBlock(A, [j + 1, j], [m - 1, m - 1]);
            C = numeric.tensor(v, numeric.dot(v, B));
            for (i = j + 1; i < m; i++) {
                Ai = A[i];
                Ci = C[i - j - 1];
                for (k = j; k < m; k++) Ai[k] -= 2 * Ci[k - j]
            }
            B = numeric.getBlock(A, [0, j + 1], [m - 1, m - 1]);
            C = numeric.tensor(numeric.dot(B, v), v);
            for (i = 0; i < m; i++) {
                Ai = A[i];
                Ci = C[i];
                for (k = j + 1; k < m; k++) Ai[k] -= 2 * Ci[k - j - 1]
            }
            B = Array(m - j - 1);
            for (i = j + 1; i < m; i++) B[i - j - 1] = Q[i];
            C = numeric.tensor(v, numeric.dot(v, B));
            for (i = j + 1; i < m; i++) {
                Qi = Q[i];
                Ci = C[i - j - 1];
                for (k = 0; k < m; k++) Qi[k] -= 2 * Ci[k]
            }
        }
    }
    return {
        H: A,
        Q: Q
    }
};
numeric.epsilon = 2.220446049250313e-16;
numeric.QRFrancis = function(H, maxiter) {
    if (typeof maxiter === "undefined") {
        maxiter = 1e4
    }
    H = numeric.clone(H);
    var H0 = numeric.clone(H);
    var s = numeric.dim(H),
        m = s[0],
        x, v, a, b, c, d, det, tr, Hloc, Q = numeric.identity(m),
        Qi, Hi, B, C, Ci, i, j, k, iter;
    if (m < 3) {
        return {
            Q: Q,
            B: [
                [0, m - 1]
            ]
        }
    }
    var epsilon = numeric.epsilon;
    for (iter = 0; iter < maxiter; iter++) {
        for (j = 0; j < m - 1; j++) {
            if (Math.abs(H[j + 1][j]) < epsilon * (Math.abs(H[j][j]) + Math.abs(H[j + 1][j + 1]))) {
                var QH1 = numeric.QRFrancis(numeric.getBlock(H, [0, 0], [j, j]), maxiter);
                var QH2 = numeric.QRFrancis(numeric.getBlock(H, [j + 1, j + 1], [m - 1, m - 1]), maxiter);
                B = Array(j + 1);
                for (i = 0; i <= j; i++) {
                    B[i] = Q[i]
                }
                C = numeric.dot(QH1.Q, B);
                for (i = 0; i <= j; i++) {
                    Q[i] = C[i]
                }
                B = Array(m - j - 1);
                for (i = j + 1; i < m; i++) {
                    B[i - j - 1] = Q[i]
                }
                C = numeric.dot(QH2.Q, B);
                for (i = j + 1; i < m; i++) {
                    Q[i] = C[i - j - 1]
                }
                return {
                    Q: Q,
                    B: QH1.B.concat(numeric.add(QH2.B, j + 1))
                }
            }
        }
        a = H[m - 2][m - 2];
        b = H[m - 2][m - 1];
        c = H[m - 1][m - 2];
        d = H[m - 1][m - 1];
        tr = a + d;
        det = a * d - b * c;
        Hloc = numeric.getBlock(H, [0, 0], [2, 2]);
        if (tr * tr >= 4 * det) {
            var s1, s2;
            s1 = .5 * (tr + Math.sqrt(tr * tr - 4 * det));
            s2 = .5 * (tr - Math.sqrt(tr * tr - 4 * det));
            Hloc = numeric.add(numeric.sub(numeric.dot(Hloc, Hloc), numeric.mul(Hloc, s1 + s2)), numeric.diag(numeric.rep([3], s1 * s2)))
        } else {
            Hloc = numeric.add(numeric.sub(numeric.dot(Hloc, Hloc), numeric.mul(Hloc, tr)), numeric.diag(numeric.rep([3], det)))
        }
        x = [Hloc[0][0], Hloc[1][0], Hloc[2][0]];
        v = numeric.house(x);
        B = [H[0], H[1], H[2]];
        C = numeric.tensor(v, numeric.dot(v, B));
        for (i = 0; i < 3; i++) {
            Hi = H[i];
            Ci = C[i];
            for (k = 0; k < m; k++) Hi[k] -= 2 * Ci[k]
        }
        B = numeric.getBlock(H, [0, 0], [m - 1, 2]);
        C = numeric.tensor(numeric.dot(B, v), v);
        for (i = 0; i < m; i++) {
            Hi = H[i];
            Ci = C[i];
            for (k = 0; k < 3; k++) Hi[k] -= 2 * Ci[k]
        }
        B = [Q[0], Q[1], Q[2]];
        C = numeric.tensor(v, numeric.dot(v, B));
        for (i = 0; i < 3; i++) {
            Qi = Q[i];
            Ci = C[i];
            for (k = 0; k < m; k++) Qi[k] -= 2 * Ci[k]
        }
        var J;
        for (j = 0; j < m - 2; j++) {
            for (k = j; k <= j + 1; k++) {
                if (Math.abs(H[k + 1][k]) < epsilon * (Math.abs(H[k][k]) + Math.abs(H[k + 1][k + 1]))) {
                    var QH1 = numeric.QRFrancis(numeric.getBlock(H, [0, 0], [k, k]), maxiter);
                    var QH2 = numeric.QRFrancis(numeric.getBlock(H, [k + 1, k + 1], [m - 1, m - 1]), maxiter);
                    B = Array(k + 1);
                    for (i = 0; i <= k; i++) {
                        B[i] = Q[i]
                    }
                    C = numeric.dot(QH1.Q, B);
                    for (i = 0; i <= k; i++) {
                        Q[i] = C[i]
                    }
                    B = Array(m - k - 1);
                    for (i = k + 1; i < m; i++) {
                        B[i - k - 1] = Q[i]
                    }
                    C = numeric.dot(QH2.Q, B);
                    for (i = k + 1; i < m; i++) {
                        Q[i] = C[i - k - 1]
                    }
                    return {
                        Q: Q,
                        B: QH1.B.concat(numeric.add(QH2.B, k + 1))
                    }
                }
            }
            J = Math.min(m - 1, j + 3);
            x = Array(J - j);
            for (i = j + 1; i <= J; i++) {
                x[i - j - 1] = H[i][j]
            }
            v = numeric.house(x);
            B = numeric.getBlock(H, [j + 1, j], [J, m - 1]);
            C = numeric.tensor(v, numeric.dot(v, B));
            for (i = j + 1; i <= J; i++) {
                Hi = H[i];
                Ci = C[i - j - 1];
                for (k = j; k < m; k++) Hi[k] -= 2 * Ci[k - j]
            }
            B = numeric.getBlock(H, [0, j + 1], [m - 1, J]);
            C = numeric.tensor(numeric.dot(B, v), v);
            for (i = 0; i < m; i++) {
                Hi = H[i];
                Ci = C[i];
                for (k = j + 1; k <= J; k++) Hi[k] -= 2 * Ci[k - j - 1]
            }
            B = Array(J - j);
            for (i = j + 1; i <= J; i++) B[i - j - 1] = Q[i];
            C = numeric.tensor(v, numeric.dot(v, B));
            for (i = j + 1; i <= J; i++) {
                Qi = Q[i];
                Ci = C[i - j - 1];
                for (k = 0; k < m; k++) Qi[k] -= 2 * Ci[k]
            }
        }
    }
    throw new Error("numeric: eigenvalue iteration does not converge -- increase maxiter?")
};
numeric.eig = function eig(A, maxiter) {
    var QH = numeric.toUpperHessenberg(A);
    var QB = numeric.QRFrancis(QH.H, maxiter);
    var T = numeric.T;
    var n = A.length,
        i, k, flag = false,
        B = QB.B,
        H = numeric.dot(QB.Q, numeric.dot(QH.H, numeric.transpose(QB.Q)));
    var Q = new T(numeric.dot(QB.Q, QH.Q)),
        Q0;
    var m = B.length,
        j;
    var a, b, c, d, p1, p2, disc, x, y, p, q, n1, n2;
    var sqrt = Math.sqrt;
    for (k = 0; k < m; k++) {
        i = B[k][0];
        if (i === B[k][1]) {} else {
            j = i + 1;
            a = H[i][i];
            b = H[i][j];
            c = H[j][i];
            d = H[j][j];
            if (b === 0 && c === 0) continue;
            p1 = -a - d;
            p2 = a * d - b * c;
            disc = p1 * p1 - 4 * p2;
            if (disc >= 0) {
                if (p1 < 0) x = -.5 * (p1 - sqrt(disc));
                else x = -.5 * (p1 + sqrt(disc));
                n1 = (a - x) * (a - x) + b * b;
                n2 = c * c + (d - x) * (d - x);
                if (n1 > n2) {
                    n1 = sqrt(n1);
                    p = (a - x) / n1;
                    q = b / n1
                } else {
                    n2 = sqrt(n2);
                    p = c / n2;
                    q = (d - x) / n2
                }
                Q0 = new T([
                    [q, -p],
                    [p, q]
                ]);
                Q.setRows(i, j, Q0.dot(Q.getRows(i, j)))
            } else {
                x = -.5 * p1;
                y = .5 * sqrt(-disc);
                n1 = (a - x) * (a - x) + b * b;
                n2 = c * c + (d - x) * (d - x);
                if (n1 > n2) {
                    n1 = sqrt(n1 + y * y);
                    p = (a - x) / n1;
                    q = b / n1;
                    x = 0;
                    y /= n1
                } else {
                    n2 = sqrt(n2 + y * y);
                    p = c / n2;
                    q = (d - x) / n2;
                    x = y / n2;
                    y = 0
                }
                Q0 = new T([
                    [q, -p],
                    [p, q]
                ], [
                    [x, y],
                    [y, -x]
                ]);
                Q.setRows(i, j, Q0.dot(Q.getRows(i, j)))
            }
        }
    }
    var R = Q.dot(A).dot(Q.transjugate()),
        n = A.length,
        E = numeric.T.identity(n);
    for (j = 0; j < n; j++) {
        if (j > 0) {
            for (k = j - 1; k >= 0; k--) {
                var Rk = R.get([k, k]),
                    Rj = R.get([j, j]);
                if (numeric.neq(Rk.x, Rj.x) || numeric.neq(Rk.y, Rj.y)) {
                    x = R.getRow(k).getBlock([k], [j - 1]);
                    y = E.getRow(j).getBlock([k], [j - 1]);
                    E.set([j, k], R.get([k, j]).neg().sub(x.dot(y)).div(Rk.sub(Rj)))
                } else {
                    E.setRow(j, E.getRow(k));
                    continue
                }
            }
        }
    }
    for (j = 0; j < n; j++) {
        x = E.getRow(j);
        E.setRow(j, x.div(x.norm2()))
    }
    E = E.transpose();
    E = Q.transjugate().dot(E);
    return {
        lambda: R.getDiag(),
        E: E
    }
};
numeric.ccsSparse = function ccsSparse(A) {
    var m = A.length,
        n, foo, i, j, counts = [];
    for (i = m - 1; i !== -1; --i) {
        foo = A[i];
        for (j in foo) {
            j = parseInt(j);
            while (j >= counts.length) counts[counts.length] = 0;
            if (foo[j] !== 0) counts[j]++
        }
    }
    var n = counts.length;
    var Ai = Array(n + 1);
    Ai[0] = 0;
    for (i = 0; i < n; ++i) Ai[i + 1] = Ai[i] + counts[i];
    var Aj = Array(Ai[n]),
        Av = Array(Ai[n]);
    for (i = m - 1; i !== -1; --i) {
        foo = A[i];
        for (j in foo) {
            if (foo[j] !== 0) {
                counts[j]--;
                Aj[Ai[j] + counts[j]] = i;
                Av[Ai[j] + counts[j]] = foo[j]
            }
        }
    }
    return [Ai, Aj, Av]
};
numeric.ccsFull = function ccsFull(A) {
    var Ai = A[0],
        Aj = A[1],
        Av = A[2],
        s = numeric.ccsDim(A),
        m = s[0],
        n = s[1],
        i, j, j0, j1, k;
    var B = numeric.rep([m, n], 0);
    for (i = 0; i < n; i++) {
        j0 = Ai[i];
        j1 = Ai[i + 1];
        for (j = j0; j < j1; ++j) {
            B[Aj[j]][i] = Av[j]
        }
    }
    return B
};
numeric.ccsTSolve = function ccsTSolve(A, b, x, bj, xj) {
    var Ai = A[0],
        Aj = A[1],
        Av = A[2],
        m = Ai.length - 1,
        max = Math.max,
        n = 0;
    if (typeof bj === "undefined") x = numeric.rep([m], 0);
    if (typeof bj === "undefined") bj = numeric.linspace(0, x.length - 1);
    if (typeof xj === "undefined") xj = [];

    function dfs(j) {
        var k;
        if (x[j] !== 0) return;
        x[j] = 1;
        for (k = Ai[j]; k < Ai[j + 1]; ++k) dfs(Aj[k]);
        xj[n] = j;
        ++n
    }
    var i, j, j0, j1, k, l, l0, l1, a;
    for (i = bj.length - 1; i !== -1; --i) {
        dfs(bj[i])
    }
    xj.length = n;
    for (i = xj.length - 1; i !== -1; --i) {
        x[xj[i]] = 0
    }
    for (i = bj.length - 1; i !== -1; --i) {
        j = bj[i];
        x[j] = b[j]
    }
    for (i = xj.length - 1; i !== -1; --i) {
        j = xj[i];
        j0 = Ai[j];
        j1 = max(Ai[j + 1], j0);
        for (k = j0; k !== j1; ++k) {
            if (Aj[k] === j) {
                x[j] /= Av[k];
                break
            }
        }
        a = x[j];
        for (k = j0; k !== j1; ++k) {
            l = Aj[k];
            if (l !== j) x[l] -= a * Av[k]
        }
    }
    return x
};
numeric.ccsDFS = function ccsDFS(n) {
    this.k = Array(n);
    this.k1 = Array(n);
    this.j = Array(n)
};
numeric.ccsDFS.prototype.dfs = function dfs(J, Ai, Aj, x, xj, Pinv) {
    var m = 0,
        foo, n = xj.length;
    var k = this.k,
        k1 = this.k1,
        j = this.j,
        km, k11;
    if (x[J] !== 0) return;
    x[J] = 1;
    j[0] = J;
    k[0] = km = Ai[J];
    k1[0] = k11 = Ai[J + 1];
    while (1) {
        if (km >= k11) {
            xj[n] = j[m];
            if (m === 0) return;
            ++n;
            --m;
            km = k[m];
            k11 = k1[m]
        } else {
            foo = Pinv[Aj[km]];
            if (x[foo] === 0) {
                x[foo] = 1;
                k[m] = km;
                ++m;
                j[m] = foo;
                km = Ai[foo];
                k1[m] = k11 = Ai[foo + 1]
            } else ++km
        }
    }
};
numeric.ccsLPSolve = function ccsLPSolve(A, B, x, xj, I, Pinv, dfs) {
    var Ai = A[0],
        Aj = A[1],
        Av = A[2],
        m = Ai.length - 1,
        n = 0;
    var Bi = B[0],
        Bj = B[1],
        Bv = B[2];
    var i, i0, i1, j, J, j0, j1, k, l, l0, l1, a;
    i0 = Bi[I];
    i1 = Bi[I + 1];
    xj.length = 0;
    for (i = i0; i < i1; ++i) {
        dfs.dfs(Pinv[Bj[i]], Ai, Aj, x, xj, Pinv)
    }
    for (i = xj.length - 1; i !== -1; --i) {
        x[xj[i]] = 0
    }
    for (i = i0; i !== i1; ++i) {
        j = Pinv[Bj[i]];
        x[j] = Bv[i]
    }
    for (i = xj.length - 1; i !== -1; --i) {
        j = xj[i];
        j0 = Ai[j];
        j1 = Ai[j + 1];
        for (k = j0; k < j1; ++k) {
            if (Pinv[Aj[k]] === j) {
                x[j] /= Av[k];
                break
            }
        }
        a = x[j];
        for (k = j0; k < j1; ++k) {
            l = Pinv[Aj[k]];
            if (l !== j) x[l] -= a * Av[k]
        }
    }
    return x
};
numeric.ccsLUP1 = function ccsLUP1(A, threshold) {
    var m = A[0].length - 1;
    var L = [numeric.rep([m + 1], 0), [],
            []
        ],
        U = [numeric.rep([m + 1], 0), [],
            []
        ];
    var Li = L[0],
        Lj = L[1],
        Lv = L[2],
        Ui = U[0],
        Uj = U[1],
        Uv = U[2];
    var x = numeric.rep([m], 0),
        xj = numeric.rep([m], 0);
    var i, j, k, j0, j1, a, e, c, d, K;
    var sol = numeric.ccsLPSolve,
        max = Math.max,
        abs = Math.abs;
    var P = numeric.linspace(0, m - 1),
        Pinv = numeric.linspace(0, m - 1);
    var dfs = new numeric.ccsDFS(m);
    if (typeof threshold === "undefined") {
        threshold = 1
    }
    for (i = 0; i < m; ++i) {
        sol(L, A, x, xj, i, Pinv, dfs);
        a = -1;
        e = -1;
        for (j = xj.length - 1; j !== -1; --j) {
            k = xj[j];
            if (k <= i) continue;
            c = abs(x[k]);
            if (c > a) {
                e = k;
                a = c
            }
        }
        if (abs(x[i]) < threshold * a) {
            j = P[i];
            a = P[e];
            P[i] = a;
            Pinv[a] = i;
            P[e] = j;
            Pinv[j] = e;
            a = x[i];
            x[i] = x[e];
            x[e] = a
        }
        a = Li[i];
        e = Ui[i];
        d = x[i];
        Lj[a] = P[i];
        Lv[a] = 1;
        ++a;
        for (j = xj.length - 1; j !== -1; --j) {
            k = xj[j];
            c = x[k];
            xj[j] = 0;
            x[k] = 0;
            if (k <= i) {
                Uj[e] = k;
                Uv[e] = c;
                ++e
            } else {
                Lj[a] = P[k];
                Lv[a] = c / d;
                ++a
            }
        }
        Li[i + 1] = a;
        Ui[i + 1] = e
    }
    for (j = Lj.length - 1; j !== -1; --j) {
        Lj[j] = Pinv[Lj[j]]
    }
    return {
        L: L,
        U: U,
        P: P,
        Pinv: Pinv
    }
};
numeric.ccsDFS0 = function ccsDFS0(n) {
    this.k = Array(n);
    this.k1 = Array(n);
    this.j = Array(n)
};
numeric.ccsDFS0.prototype.dfs = function dfs(J, Ai, Aj, x, xj, Pinv, P) {
    var m = 0,
        foo, n = xj.length;
    var k = this.k,
        k1 = this.k1,
        j = this.j,
        km, k11;
    if (x[J] !== 0) return;
    x[J] = 1;
    j[0] = J;
    k[0] = km = Ai[Pinv[J]];
    k1[0] = k11 = Ai[Pinv[J] + 1];
    while (1) {
        if (isNaN(km)) throw new Error("Ow!");
        if (km >= k11) {
            xj[n] = Pinv[j[m]];
            if (m === 0) return;
            ++n;
            --m;
            km = k[m];
            k11 = k1[m]
        } else {
            foo = Aj[km];
            if (x[foo] === 0) {
                x[foo] = 1;
                k[m] = km;
                ++m;
                j[m] = foo;
                foo = Pinv[foo];
                km = Ai[foo];
                k1[m] = k11 = Ai[foo + 1]
            } else ++km
        }
    }
};
numeric.ccsLPSolve0 = function ccsLPSolve0(A, B, y, xj, I, Pinv, P, dfs) {
    var Ai = A[0],
        Aj = A[1],
        Av = A[2],
        m = Ai.length - 1,
        n = 0;
    var Bi = B[0],
        Bj = B[1],
        Bv = B[2];
    var i, i0, i1, j, J, j0, j1, k, l, l0, l1, a;
    i0 = Bi[I];
    i1 = Bi[I + 1];
    xj.length = 0;
    for (i = i0; i < i1; ++i) {
        dfs.dfs(Bj[i], Ai, Aj, y, xj, Pinv, P)
    }
    for (i = xj.length - 1; i !== -1; --i) {
        j = xj[i];
        y[P[j]] = 0
    }
    for (i = i0; i !== i1; ++i) {
        j = Bj[i];
        y[j] = Bv[i]
    }
    for (i = xj.length - 1; i !== -1; --i) {
        j = xj[i];
        l = P[j];
        j0 = Ai[j];
        j1 = Ai[j + 1];
        for (k = j0; k < j1; ++k) {
            if (Aj[k] === l) {
                y[l] /= Av[k];
                break
            }
        }
        a = y[l];
        for (k = j0; k < j1; ++k) y[Aj[k]] -= a * Av[k];
        y[l] = a
    }
};
numeric.ccsLUP0 = function ccsLUP0(A, threshold) {
    var m = A[0].length - 1;
    var L = [numeric.rep([m + 1], 0), [],
            []
        ],
        U = [numeric.rep([m + 1], 0), [],
            []
        ];
    var Li = L[0],
        Lj = L[1],
        Lv = L[2],
        Ui = U[0],
        Uj = U[1],
        Uv = U[2];
    var y = numeric.rep([m], 0),
        xj = numeric.rep([m], 0);
    var i, j, k, j0, j1, a, e, c, d, K;
    var sol = numeric.ccsLPSolve0,
        max = Math.max,
        abs = Math.abs;
    var P = numeric.linspace(0, m - 1),
        Pinv = numeric.linspace(0, m - 1);
    var dfs = new numeric.ccsDFS0(m);
    if (typeof threshold === "undefined") {
        threshold = 1
    }
    for (i = 0; i < m; ++i) {
        sol(L, A, y, xj, i, Pinv, P, dfs);
        a = -1;
        e = -1;
        for (j = xj.length - 1; j !== -1; --j) {
            k = xj[j];
            if (k <= i) continue;
            c = abs(y[P[k]]);
            if (c > a) {
                e = k;
                a = c
            }
        }
        if (abs(y[P[i]]) < threshold * a) {
            j = P[i];
            a = P[e];
            P[i] = a;
            Pinv[a] = i;
            P[e] = j;
            Pinv[j] = e
        }
        a = Li[i];
        e = Ui[i];
        d = y[P[i]];
        Lj[a] = P[i];
        Lv[a] = 1;
        ++a;
        for (j = xj.length - 1; j !== -1; --j) {
            k = xj[j];
            c = y[P[k]];
            xj[j] = 0;
            y[P[k]] = 0;
            if (k <= i) {
                Uj[e] = k;
                Uv[e] = c;
                ++e
            } else {
                Lj[a] = P[k];
                Lv[a] = c / d;
                ++a
            }
        }
        Li[i + 1] = a;
        Ui[i + 1] = e
    }
    for (j = Lj.length - 1; j !== -1; --j) {
        Lj[j] = Pinv[Lj[j]]
    }
    return {
        L: L,
        U: U,
        P: P,
        Pinv: Pinv
    }
};
numeric.ccsLUP = numeric.ccsLUP0;
numeric.ccsDim = function ccsDim(A) {
    return [numeric.sup(A[1]) + 1, A[0].length - 1]
};
numeric.ccsGetBlock = function ccsGetBlock(A, i, j) {
    var s = numeric.ccsDim(A),
        m = s[0],
        n = s[1];
    if (typeof i === "undefined") {
        i = numeric.linspace(0, m - 1)
    } else if (typeof i === "number") {
        i = [i]
    }
    if (typeof j === "undefined") {
        j = numeric.linspace(0, n - 1)
    } else if (typeof j === "number") {
        j = [j]
    }
    var p, p0, p1, P = i.length,
        q, Q = j.length,
        r, jq, ip;
    var Bi = numeric.rep([n], 0),
        Bj = [],
        Bv = [],
        B = [Bi, Bj, Bv];
    var Ai = A[0],
        Aj = A[1],
        Av = A[2];
    var x = numeric.rep([m], 0),
        count = 0,
        flags = numeric.rep([m], 0);
    for (q = 0; q < Q; ++q) {
        jq = j[q];
        var q0 = Ai[jq];
        var q1 = Ai[jq + 1];
        for (p = q0; p < q1; ++p) {
            r = Aj[p];
            flags[r] = 1;
            x[r] = Av[p]
        }
        for (p = 0; p < P; ++p) {
            ip = i[p];
            if (flags[ip]) {
                Bj[count] = p;
                Bv[count] = x[i[p]];
                ++count
            }
        }
        for (p = q0; p < q1; ++p) {
            r = Aj[p];
            flags[r] = 0
        }
        Bi[q + 1] = count
    }
    return B
};
numeric.ccsDot = function ccsDot(A, B) {
    var Ai = A[0],
        Aj = A[1],
        Av = A[2];
    var Bi = B[0],
        Bj = B[1],
        Bv = B[2];
    var sA = numeric.ccsDim(A),
        sB = numeric.ccsDim(B);
    var m = sA[0],
        n = sA[1],
        o = sB[1];
    var x = numeric.rep([m], 0),
        flags = numeric.rep([m], 0),
        xj = Array(m);
    var Ci = numeric.rep([o], 0),
        Cj = [],
        Cv = [],
        C = [Ci, Cj, Cv];
    var i, j, k, j0, j1, i0, i1, l, p, a, b;
    for (k = 0; k !== o; ++k) {
        j0 = Bi[k];
        j1 = Bi[k + 1];
        p = 0;
        for (j = j0; j < j1; ++j) {
            a = Bj[j];
            b = Bv[j];
            i0 = Ai[a];
            i1 = Ai[a + 1];
            for (i = i0; i < i1; ++i) {
                l = Aj[i];
                if (flags[l] === 0) {
                    xj[p] = l;
                    flags[l] = 1;
                    p = p + 1
                }
                x[l] = x[l] + Av[i] * b
            }
        }
        j0 = Ci[k];
        j1 = j0 + p;
        Ci[k + 1] = j1;
        for (j = p - 1; j !== -1; --j) {
            b = j0 + j;
            i = xj[j];
            Cj[b] = i;
            Cv[b] = x[i];
            flags[i] = 0;
            x[i] = 0
        }
        Ci[k + 1] = Ci[k] + p
    }
    return C
};
numeric.ccsLUPSolve = function ccsLUPSolve(LUP, B) {
    var L = LUP.L,
        U = LUP.U,
        P = LUP.P;
    var Bi = B[0];
    var flag = false;
    if (typeof Bi !== "object") {
        B = [
            [0, B.length], numeric.linspace(0, B.length - 1), B
        ];
        Bi = B[0];
        flag = true
    }
    var Bj = B[1],
        Bv = B[2];
    var n = L[0].length - 1,
        m = Bi.length - 1;
    var x = numeric.rep([n], 0),
        xj = Array(n);
    var b = numeric.rep([n], 0),
        bj = Array(n);
    var Xi = numeric.rep([m + 1], 0),
        Xj = [],
        Xv = [];
    var sol = numeric.ccsTSolve;
    var i, j, j0, j1, k, J, N = 0;
    for (i = 0; i < m; ++i) {
        k = 0;
        j0 = Bi[i];
        j1 = Bi[i + 1];
        for (j = j0; j < j1; ++j) {
            J = LUP.Pinv[Bj[j]];
            bj[k] = J;
            b[J] = Bv[j];
            ++k
        }
        bj.length = k;
        sol(L, b, x, bj, xj);
        for (j = bj.length - 1; j !== -1; --j) b[bj[j]] = 0;
        sol(U, x, b, xj, bj);
        if (flag) return b;
        for (j = xj.length - 1; j !== -1; --j) x[xj[j]] = 0;
        for (j = bj.length - 1; j !== -1; --j) {
            J = bj[j];
            Xj[N] = J;
            Xv[N] = b[J];
            b[J] = 0;
            ++N
        }
        Xi[i + 1] = N
    }
    return [Xi, Xj, Xv]
};
numeric.ccsbinop = function ccsbinop(body, setup) {
    if (typeof setup === "undefined") setup = "";
    return Function("X", "Y", "var Xi = X[0], Xj = X[1], Xv = X[2];\n" + "var Yi = Y[0], Yj = Y[1], Yv = Y[2];\n" + "var n = Xi.length-1,m = Math.max(numeric.sup(Xj),numeric.sup(Yj))+1;\n" + "var Zi = numeric.rep([n+1],0), Zj = [], Zv = [];\n" + "var x = numeric.rep([m],0),y = numeric.rep([m],0);\n" + "var xk,yk,zk;\n" + "var i,j,j0,j1,k,p=0;\n" + setup + "for(i=0;i<n;++i) {\n" + "  j0 = Xi[i]; j1 = Xi[i+1];\n" + "  for(j=j0;j!==j1;++j) {\n" + "    k = Xj[j];\n" + "    x[k] = 1;\n" + "    Zj[p] = k;\n" + "    ++p;\n" + "  }\n" + "  j0 = Yi[i]; j1 = Yi[i+1];\n" + "  for(j=j0;j!==j1;++j) {\n" + "    k = Yj[j];\n" + "    y[k] = Yv[j];\n" + "    if(x[k] === 0) {\n" + "      Zj[p] = k;\n" + "      ++p;\n" + "    }\n" + "  }\n" + "  Zi[i+1] = p;\n" + "  j0 = Xi[i]; j1 = Xi[i+1];\n" + "  for(j=j0;j!==j1;++j) x[Xj[j]] = Xv[j];\n" + "  j0 = Zi[i]; j1 = Zi[i+1];\n" + "  for(j=j0;j!==j1;++j) {\n" + "    k = Zj[j];\n" + "    xk = x[k];\n" + "    yk = y[k];\n" + body + "\n" + "    Zv[j] = zk;\n" + "  }\n" + "  j0 = Xi[i]; j1 = Xi[i+1];\n" + "  for(j=j0;j!==j1;++j) x[Xj[j]] = 0;\n" + "  j0 = Yi[i]; j1 = Yi[i+1];\n" + "  for(j=j0;j!==j1;++j) y[Yj[j]] = 0;\n" + "}\n" + "return [Zi,Zj,Zv];")
};
(function() {
    var k, A, B, C;
    for (k in numeric.ops2) {
        if (isFinite(eval("1" + numeric.ops2[k] + "0"))) A = "[Y[0],Y[1],numeric." + k + "(X,Y[2])]";
        else A = "NaN";
        if (isFinite(eval("0" + numeric.ops2[k] + "1"))) B = "[X[0],X[1],numeric." + k + "(X[2],Y)]";
        else B = "NaN";
        if (isFinite(eval("1" + numeric.ops2[k] + "0")) && isFinite(eval("0" + numeric.ops2[k] + "1"))) C = "numeric.ccs" + k + "MM(X,Y)";
        else C = "NaN";
        numeric["ccs" + k + "MM"] = numeric.ccsbinop("zk = xk " + numeric.ops2[k] + "yk;");
        numeric["ccs" + k] = Function("X", "Y", 'if(typeof X === "number") return ' + A + ";\n" + 'if(typeof Y === "number") return ' + B + ";\n" + "return " + C + ";\n")
    }
})();
numeric.ccsScatter = function ccsScatter(A) {
    var Ai = A[0],
        Aj = A[1],
        Av = A[2];
    var n = numeric.sup(Aj) + 1,
        m = Ai.length;
    var Ri = numeric.rep([n], 0),
        Rj = Array(m),
        Rv = Array(m);
    var counts = numeric.rep([n], 0),
        i;
    for (i = 0; i < m; ++i) counts[Aj[i]]++;
    for (i = 0; i < n; ++i) Ri[i + 1] = Ri[i] + counts[i];
    var ptr = Ri.slice(0),
        k, Aii;
    for (i = 0; i < m; ++i) {
        Aii = Aj[i];
        k = ptr[Aii];
        Rj[k] = Ai[i];
        Rv[k] = Av[i];
        ptr[Aii] = ptr[Aii] + 1
    }
    return [Ri, Rj, Rv]
};
numeric.ccsGather = function ccsGather(A) {
    var Ai = A[0],
        Aj = A[1],
        Av = A[2];
    var n = Ai.length - 1,
        m = Aj.length;
    var Ri = Array(m),
        Rj = Array(m),
        Rv = Array(m);
    var i, j, j0, j1, p;
    p = 0;
    for (i = 0; i < n; ++i) {
        j0 = Ai[i];
        j1 = Ai[i + 1];
        for (j = j0; j !== j1; ++j) {
            Rj[p] = i;
            Ri[p] = Aj[j];
            Rv[p] = Av[j];
            ++p
        }
    }
    return [Ri, Rj, Rv]
};
numeric.sdim = function dim(A, ret, k) {
    if (typeof ret === "undefined") {
        ret = []
    }
    if (typeof A !== "object") return ret;
    if (typeof k === "undefined") {
        k = 0
    }
    if (!(k in ret)) {
        ret[k] = 0
    }
    if (A.length > ret[k]) ret[k] = A.length;
    var i;
    for (i in A) {
        if (A.hasOwnProperty(i)) dim(A[i], ret, k + 1)
    }
    return ret
};
numeric.sclone = function clone(A, k, n) {
    if (typeof k === "undefined") {
        k = 0
    }
    if (typeof n === "undefined") {
        n = numeric.sdim(A).length
    }
    var i, ret = Array(A.length);
    if (k === n - 1) {
        for (i in A) {
            if (A.hasOwnProperty(i)) ret[i] = A[i]
        }
        return ret
    }
    for (i in A) {
        if (A.hasOwnProperty(i)) ret[i] = clone(A[i], k + 1, n)
    }
    return ret
};
numeric.sdiag = function diag(d) {
    var n = d.length,
        i, ret = Array(n),
        i1, i2, i3;
    for (i = n - 1; i >= 1; i -= 2) {
        i1 = i - 1;
        ret[i] = [];
        ret[i][i] = d[i];
        ret[i1] = [];
        ret[i1][i1] = d[i1]
    }
    if (i === 0) {
        ret[0] = [];
        ret[0][0] = d[i]
    }
    return ret
};
numeric.sidentity = function identity(n) {
    return numeric.sdiag(numeric.rep([n], 1))
};
numeric.stranspose = function transpose(A) {
    var ret = [],
        n = A.length,
        i, j, Ai;
    for (i in A) {
        if (!A.hasOwnProperty(i)) continue;
        Ai = A[i];
        for (j in Ai) {
            if (!Ai.hasOwnProperty(j)) continue;
            if (typeof ret[j] !== "object") {
                ret[j] = []
            }
            ret[j][i] = Ai[j]
        }
    }
    return ret
};
numeric.sLUP = function LUP(A, tol) {
    throw new Error("The function numeric.sLUP had a bug in it and has been removed. Please use the new numeric.ccsLUP function instead.")
};
numeric.sdotMM = function dotMM(A, B) {
    var p = A.length,
        q = B.length,
        BT = numeric.stranspose(B),
        r = BT.length,
        Ai, BTk;
    var i, j, k, accum;
    var ret = Array(p),
        reti;
    for (i = p - 1; i >= 0; i--) {
        reti = [];
        Ai = A[i];
        for (k = r - 1; k >= 0; k--) {
            accum = 0;
            BTk = BT[k];
            for (j in Ai) {
                if (!Ai.hasOwnProperty(j)) continue;
                if (j in BTk) {
                    accum += Ai[j] * BTk[j]
                }
            }
            if (accum) reti[k] = accum
        }
        ret[i] = reti
    }
    return ret
};
numeric.sdotMV = function dotMV(A, x) {
    var p = A.length,
        Ai, i, j;
    var ret = Array(p),
        accum;
    for (i = p - 1; i >= 0; i--) {
        Ai = A[i];
        accum = 0;
        for (j in Ai) {
            if (!Ai.hasOwnProperty(j)) continue;
            if (x[j]) accum += Ai[j] * x[j]
        }
        if (accum) ret[i] = accum
    }
    return ret
};
numeric.sdotVM = function dotMV(x, A) {
    var i, j, Ai, alpha;
    var ret = [],
        accum;
    for (i in x) {
        if (!x.hasOwnProperty(i)) continue;
        Ai = A[i];
        alpha = x[i];
        for (j in Ai) {
            if (!Ai.hasOwnProperty(j)) continue;
            if (!ret[j]) {
                ret[j] = 0
            }
            ret[j] += alpha * Ai[j]
        }
    }
    return ret
};
numeric.sdotVV = function dotVV(x, y) {
    var i, ret = 0;
    for (i in x) {
        if (x[i] && y[i]) ret += x[i] * y[i]
    }
    return ret
};
numeric.sdot = function dot(A, B) {
    var m = numeric.sdim(A).length,
        n = numeric.sdim(B).length;
    var k = m * 1e3 + n;
    switch (k) {
        case 0:
            return A * B;
        case 1001:
            return numeric.sdotVV(A, B);
        case 2001:
            return numeric.sdotMV(A, B);
        case 1002:
            return numeric.sdotVM(A, B);
        case 2002:
            return numeric.sdotMM(A, B);
        default:
            throw new Error("numeric.sdot not implemented for tensors of order " + m + " and " + n)
    }
};
numeric.sscatter = function scatter(V) {
    var n = V[0].length,
        Vij, i, j, m = V.length,
        A = [],
        Aj;
    for (i = n - 1; i >= 0; --i) {
        if (!V[m - 1][i]) continue;
        Aj = A;
        for (j = 0; j < m - 2; j++) {
            Vij = V[j][i];
            if (!Aj[Vij]) Aj[Vij] = [];
            Aj = Aj[Vij]
        }
        Aj[V[j][i]] = V[j + 1][i]
    }
    return A
};
numeric.sgather = function gather(A, ret, k) {
    if (typeof ret === "undefined") ret = [];
    if (typeof k === "undefined") k = [];
    var n, i, Ai;
    n = k.length;
    for (i in A) {
        if (A.hasOwnProperty(i)) {
            k[n] = parseInt(i);
            Ai = A[i];
            if (typeof Ai === "number") {
                if (Ai) {
                    if (ret.length === 0) {
                        for (i = n + 1; i >= 0; --i) ret[i] = []
                    }
                    for (i = n; i >= 0; --i) ret[i].push(k[i]);
                    ret[n + 1].push(Ai)
                }
            } else gather(Ai, ret, k)
        }
    }
    if (k.length > n) k.pop();
    return ret
};
numeric.cLU = function LU(A) {
    var I = A[0],
        J = A[1],
        V = A[2];
    var p = I.length,
        m = 0,
        i, j, k, a, b, c;
    for (i = 0; i < p; i++)
        if (I[i] > m) m = I[i];
    m++;
    var L = Array(m),
        U = Array(m),
        left = numeric.rep([m], Infinity),
        right = numeric.rep([m], -Infinity);
    var Ui, Uj, alpha;
    for (k = 0; k < p; k++) {
        i = I[k];
        j = J[k];
        if (j < left[i]) left[i] = j;
        if (j > right[i]) right[i] = j
    }
    for (i = 0; i < m - 1; i++) {
        if (right[i] > right[i + 1]) right[i + 1] = right[i]
    }
    for (i = m - 1; i >= 1; i--) {
        if (left[i] < left[i - 1]) left[i - 1] = left[i]
    }
    var countL = 0,
        countU = 0;
    for (i = 0; i < m; i++) {
        U[i] = numeric.rep([right[i] - left[i] + 1], 0);
        L[i] = numeric.rep([i - left[i]], 0);
        countL += i - left[i] + 1;
        countU += right[i] - i + 1
    }
    for (k = 0; k < p; k++) {
        i = I[k];
        U[i][J[k] - left[i]] = V[k]
    }
    for (i = 0; i < m - 1; i++) {
        a = i - left[i];
        Ui = U[i];
        for (j = i + 1; left[j] <= i && j < m; j++) {
            b = i - left[j];
            c = right[i] - i;
            Uj = U[j];
            alpha = Uj[b] / Ui[a];
            if (alpha) {
                for (k = 1; k <= c; k++) {
                    Uj[k + b] -= alpha * Ui[k + a]
                }
                L[j][i - left[j]] = alpha
            }
        }
    }
    var Ui = [],
        Uj = [],
        Uv = [],
        Li = [],
        Lj = [],
        Lv = [];
    var p, q, foo;
    p = 0;
    q = 0;
    for (i = 0; i < m; i++) {
        a = left[i];
        b = right[i];
        foo = U[i];
        for (j = i; j <= b; j++) {
            if (foo[j - a]) {
                Ui[p] = i;
                Uj[p] = j;
                Uv[p] = foo[j - a];
                p++
            }
        }
        foo = L[i];
        for (j = a; j < i; j++) {
            if (foo[j - a]) {
                Li[q] = i;
                Lj[q] = j;
                Lv[q] = foo[j - a];
                q++
            }
        }
        Li[q] = i;
        Lj[q] = i;
        Lv[q] = 1;
        q++
    }
    return {
        U: [Ui, Uj, Uv],
        L: [Li, Lj, Lv]
    }
};
numeric.cLUsolve = function LUsolve(lu, b) {
    var L = lu.L,
        U = lu.U,
        ret = numeric.clone(b);
    var Li = L[0],
        Lj = L[1],
        Lv = L[2];
    var Ui = U[0],
        Uj = U[1],
        Uv = U[2];
    var p = Ui.length,
        q = Li.length;
    var m = ret.length,
        i, j, k;
    k = 0;
    for (i = 0; i < m; i++) {
        while (Lj[k] < i) {
            ret[i] -= Lv[k] * ret[Lj[k]];
            k++
        }
        k++
    }
    k = p - 1;
    for (i = m - 1; i >= 0; i--) {
        while (Uj[k] > i) {
            ret[i] -= Uv[k] * ret[Uj[k]];
            k--
        }
        ret[i] /= Uv[k];
        k--
    }
    return ret
};
numeric.cgrid = function grid(n, shape) {
    if (typeof n === "number") n = [n, n];
    var ret = numeric.rep(n, -1);
    var i, j, count;
    if (typeof shape !== "function") {
        switch (shape) {
            case "L":
                shape = function(i, j) {
                    return i >= n[0] / 2 || j < n[1] / 2
                };
                break;
            default:
                shape = function(i, j) {
                    return true
                };
                break
        }
    }
    count = 0;
    for (i = 1; i < n[0] - 1; i++)
        for (j = 1; j < n[1] - 1; j++)
            if (shape(i, j)) {
                ret[i][j] = count;
                count++
            }
    return ret
};
numeric.cdelsq = function delsq(g) {
    var dir = [
        [-1, 0],
        [0, -1],
        [0, 1],
        [1, 0]
    ];
    var s = numeric.dim(g),
        m = s[0],
        n = s[1],
        i, j, k, p, q;
    var Li = [],
        Lj = [],
        Lv = [];
    for (i = 1; i < m - 1; i++)
        for (j = 1; j < n - 1; j++) {
            if (g[i][j] < 0) continue;
            for (k = 0; k < 4; k++) {
                p = i + dir[k][0];
                q = j + dir[k][1];
                if (g[p][q] < 0) continue;
                Li.push(g[i][j]);
                Lj.push(g[p][q]);
                Lv.push(-1)
            }
            Li.push(g[i][j]);
            Lj.push(g[i][j]);
            Lv.push(4)
        }
    return [Li, Lj, Lv]
};
numeric.cdotMV = function dotMV(A, x) {
    var ret, Ai = A[0],
        Aj = A[1],
        Av = A[2],
        k, p = Ai.length,
        N;
    N = 0;
    for (k = 0; k < p; k++) {
        if (Ai[k] > N) N = Ai[k]
    }
    N++;
    ret = numeric.rep([N], 0);
    for (k = 0; k < p; k++) {
        ret[Ai[k]] += Av[k] * x[Aj[k]]
    }
    return ret
};
numeric.Spline = function Spline(x, yl, yr, kl, kr) {
    this.x = x;
    this.yl = yl;
    this.yr = yr;
    this.kl = kl;
    this.kr = kr
};
numeric.Spline.prototype._at = function _at(x1, p) {
    var x = this.x;
    var yl = this.yl;
    var yr = this.yr;
    var kl = this.kl;
    var kr = this.kr;
    var x1, a, b, t;
    var add = numeric.add,
        sub = numeric.sub,
        mul = numeric.mul;
    a = sub(mul(kl[p], x[p + 1] - x[p]), sub(yr[p + 1], yl[p]));
    b = add(mul(kr[p + 1], x[p] - x[p + 1]), sub(yr[p + 1], yl[p]));
    t = (x1 - x[p]) / (x[p + 1] - x[p]);
    var s = t * (1 - t);
    return add(add(add(mul(1 - t, yl[p]), mul(t, yr[p + 1])), mul(a, s * (1 - t))), mul(b, s * t))
};
numeric.Spline.prototype.at = function at(x0) {
    if (typeof x0 === "number") {
        var x = this.x;
        var n = x.length;
        var p, q, mid, floor = Math.floor,
            a, b, t;
        p = 0;
        q = n - 1;
        while (q - p > 1) {
            mid = floor((p + q) / 2);
            if (x[mid] <= x0) p = mid;
            else q = mid
        }
        return this._at(x0, p)
    }
    var n = x0.length,
        i, ret = Array(n);
    for (i = n - 1; i !== -1; --i) ret[i] = this.at(x0[i]);
    return ret
};
numeric.Spline.prototype.diff = function diff() {
    var x = this.x;
    var yl = this.yl;
    var yr = this.yr;
    var kl = this.kl;
    var kr = this.kr;
    var n = yl.length;
    var i, dx, dy;
    var zl = kl,
        zr = kr,
        pl = Array(n),
        pr = Array(n);
    var add = numeric.add,
        mul = numeric.mul,
        div = numeric.div,
        sub = numeric.sub;
    for (i = n - 1; i !== -1; --i) {
        dx = x[i + 1] - x[i];
        dy = sub(yr[i + 1], yl[i]);
        pl[i] = div(add(mul(dy, 6), mul(kl[i], -4 * dx), mul(kr[i + 1], -2 * dx)), dx * dx);
        pr[i + 1] = div(add(mul(dy, -6), mul(kl[i], 2 * dx), mul(kr[i + 1], 4 * dx)), dx * dx)
    }
    return new numeric.Spline(x, zl, zr, pl, pr)
};
numeric.Spline.prototype.roots = function roots() {
    function sqr(x) {
        return x * x
    }

    function heval(y0, y1, k0, k1, x) {
        var A = k0 * 2 - (y1 - y0);
        var B = -k1 * 2 + (y1 - y0);
        var t = (x + 1) * .5;
        var s = t * (1 - t);
        return (1 - t) * y0 + t * y1 + A * s * (1 - t) + B * s * t
    }
    var ret = [];
    var x = this.x,
        yl = this.yl,
        yr = this.yr,
        kl = this.kl,
        kr = this.kr;
    if (typeof yl[0] === "number") {
        yl = [yl];
        yr = [yr];
        kl = [kl];
        kr = [kr]
    }
    var m = yl.length,
        n = x.length - 1,
        i, j, k, y, s, t;
    var ai, bi, ci, di, ret = Array(m),
        ri, k0, k1, y0, y1, A, B, D, dx, cx, stops, z0, z1, zm, t0, t1, tm;
    var sqrt = Math.sqrt;
    for (i = 0; i !== m; ++i) {
        ai = yl[i];
        bi = yr[i];
        ci = kl[i];
        di = kr[i];
        ri = [];
        for (j = 0; j !== n; j++) {
            if (j > 0 && bi[j] * ai[j] < 0) ri.push(x[j]);
            dx = x[j + 1] - x[j];
            cx = x[j];
            y0 = ai[j];
            y1 = bi[j + 1];
            k0 = ci[j] / dx;
            k1 = di[j + 1] / dx;
            D = sqr(k0 - k1 + 3 * (y0 - y1)) + 12 * k1 * y0;
            A = k1 + 3 * y0 + 2 * k0 - 3 * y1;
            B = 3 * (k1 + k0 + 2 * (y0 - y1));
            if (D <= 0) {
                z0 = A / B;
                if (z0 > x[j] && z0 < x[j + 1]) stops = [x[j], z0, x[j + 1]];
                else stops = [x[j], x[j + 1]]
            } else {
                z0 = (A - sqrt(D)) / B;
                z1 = (A + sqrt(D)) / B;
                stops = [x[j]];
                if (z0 > x[j] && z0 < x[j + 1]) stops.push(z0);
                if (z1 > x[j] && z1 < x[j + 1]) stops.push(z1);
                stops.push(x[j + 1])
            }
            t0 = stops[0];
            z0 = this._at(t0, j);
            for (k = 0; k < stops.length - 1; k++) {
                t1 = stops[k + 1];
                z1 = this._at(t1, j);
                if (z0 === 0) {
                    ri.push(t0);
                    t0 = t1;
                    z0 = z1;
                    continue
                }
                if (z1 === 0 || z0 * z1 > 0) {
                    t0 = t1;
                    z0 = z1;
                    continue
                }
                var side = 0;
                while (1) {
                    tm = (z0 * t1 - z1 * t0) / (z0 - z1);
                    if (tm <= t0 || tm >= t1) {
                        break
                    }
                    zm = this._at(tm, j);
                    if (zm * z1 > 0) {
                        t1 = tm;
                        z1 = zm;
                        if (side === -1) z0 *= .5;
                        side = -1
                    } else if (zm * z0 > 0) {
                        t0 = tm;
                        z0 = zm;
                        if (side === 1) z1 *= .5;
                        side = 1
                    } else break
                }
                ri.push(tm);
                t0 = stops[k + 1];
                z0 = this._at(t0, j)
            }
            if (z1 === 0) ri.push(t1)
        }
        ret[i] = ri
    }
    if (typeof this.yl[0] === "number") return ret[0];
    return ret
};
numeric.spline = function spline(x, y, k1, kn) {
    var n = x.length,
        b = [],
        dx = [],
        dy = [];
    var i;
    var sub = numeric.sub,
        mul = numeric.mul,
        add = numeric.add;
    for (i = n - 2; i >= 0; i--) {
        dx[i] = x[i + 1] - x[i];
        dy[i] = sub(y[i + 1], y[i])
    }
    if (typeof k1 === "string" || typeof kn === "string") {
        k1 = kn = "periodic"
    }
    var T = [
        [],
        [],
        []
    ];
    switch (typeof k1) {
        case "undefined":
            b[0] = mul(3 / (dx[0] * dx[0]), dy[0]);
            T[0].push(0, 0);
            T[1].push(0, 1);
            T[2].push(2 / dx[0], 1 / dx[0]);
            break;
        case "string":
            b[0] = add(mul(3 / (dx[n - 2] * dx[n - 2]), dy[n - 2]), mul(3 / (dx[0] * dx[0]), dy[0]));
            T[0].push(0, 0, 0);
            T[1].push(n - 2, 0, 1);
            T[2].push(1 / dx[n - 2], 2 / dx[n - 2] + 2 / dx[0], 1 / dx[0]);
            break;
        default:
            b[0] = k1;
            T[0].push(0);
            T[1].push(0);
            T[2].push(1);
            break
    }
    for (i = 1; i < n - 1; i++) {
        b[i] = add(mul(3 / (dx[i - 1] * dx[i - 1]), dy[i - 1]), mul(3 / (dx[i] * dx[i]), dy[i]));
        T[0].push(i, i, i);
        T[1].push(i - 1, i, i + 1);
        T[2].push(1 / dx[i - 1], 2 / dx[i - 1] + 2 / dx[i], 1 / dx[i])
    }
    switch (typeof kn) {
        case "undefined":
            b[n - 1] = mul(3 / (dx[n - 2] * dx[n - 2]), dy[n - 2]);
            T[0].push(n - 1, n - 1);
            T[1].push(n - 2, n - 1);
            T[2].push(1 / dx[n - 2], 2 / dx[n - 2]);
            break;
        case "string":
            T[1][T[1].length - 1] = 0;
            break;
        default:
            b[n - 1] = kn;
            T[0].push(n - 1);
            T[1].push(n - 1);
            T[2].push(1);
            break
    }
    if (typeof b[0] !== "number") b = numeric.transpose(b);
    else b = [b];
    var k = Array(b.length);
    if (typeof k1 === "string") {
        for (i = k.length - 1; i !== -1; --i) {
            k[i] = numeric.ccsLUPSolve(numeric.ccsLUP(numeric.ccsScatter(T)), b[i]);
            k[i][n - 1] = k[i][0]
        }
    } else {
        for (i = k.length - 1; i !== -1; --i) {
            k[i] = numeric.cLUsolve(numeric.cLU(T), b[i])
        }
    }
    if (typeof y[0] === "number") k = k[0];
    else k = numeric.transpose(k);
    return new numeric.Spline(x, y, y, k, k)
};
numeric.fftpow2 = function fftpow2(x, y) {
    var n = x.length;
    if (n === 1) return;
    var cos = Math.cos,
        sin = Math.sin,
        i, j;
    var xe = Array(n / 2),
        ye = Array(n / 2),
        xo = Array(n / 2),
        yo = Array(n / 2);
    j = n / 2;
    for (i = n - 1; i !== -1; --i) {
        --j;
        xo[j] = x[i];
        yo[j] = y[i];
        --i;
        xe[j] = x[i];
        ye[j] = y[i]
    }
    fftpow2(xe, ye);
    fftpow2(xo, yo);
    j = n / 2;
    var t, k = -6.283185307179586 / n,
        ci, si;
    for (i = n - 1; i !== -1; --i) {
        --j;
        if (j === -1) j = n / 2 - 1;
        t = k * i;
        ci = cos(t);
        si = sin(t);
        x[i] = xe[j] + ci * xo[j] - si * yo[j];
        y[i] = ye[j] + ci * yo[j] + si * xo[j]
    }
};
numeric._ifftpow2 = function _ifftpow2(x, y) {
    var n = x.length;
    if (n === 1) return;
    var cos = Math.cos,
        sin = Math.sin,
        i, j;
    var xe = Array(n / 2),
        ye = Array(n / 2),
        xo = Array(n / 2),
        yo = Array(n / 2);
    j = n / 2;
    for (i = n - 1; i !== -1; --i) {
        --j;
        xo[j] = x[i];
        yo[j] = y[i];
        --i;
        xe[j] = x[i];
        ye[j] = y[i]
    }
    _ifftpow2(xe, ye);
    _ifftpow2(xo, yo);
    j = n / 2;
    var t, k = 6.283185307179586 / n,
        ci, si;
    for (i = n - 1; i !== -1; --i) {
        --j;
        if (j === -1) j = n / 2 - 1;
        t = k * i;
        ci = cos(t);
        si = sin(t);
        x[i] = xe[j] + ci * xo[j] - si * yo[j];
        y[i] = ye[j] + ci * yo[j] + si * xo[j]
    }
};
numeric.ifftpow2 = function ifftpow2(x, y) {
    numeric._ifftpow2(x, y);
    numeric.diveq(x, x.length);
    numeric.diveq(y, y.length)
};
numeric.convpow2 = function convpow2(ax, ay, bx, by) {
    numeric.fftpow2(ax, ay);
    numeric.fftpow2(bx, by);
    var i, n = ax.length,
        axi, bxi, ayi, byi;
    for (i = n - 1; i !== -1; --i) {
        axi = ax[i];
        ayi = ay[i];
        bxi = bx[i];
        byi = by[i];
        ax[i] = axi * bxi - ayi * byi;
        ay[i] = axi * byi + ayi * bxi
    }
    numeric.ifftpow2(ax, ay)
};
numeric.T.prototype.fft = function fft() {
    var x = this.x,
        y = this.y;
    var n = x.length,
        log = Math.log,
        log2 = log(2),
        p = Math.ceil(log(2 * n - 1) / log2),
        m = Math.pow(2, p);
    var cx = numeric.rep([m], 0),
        cy = numeric.rep([m], 0),
        cos = Math.cos,
        sin = Math.sin;
    var k, c = -3.141592653589793 / n,
        t;
    var a = numeric.rep([m], 0),
        b = numeric.rep([m], 0),
        nhalf = Math.floor(n / 2);
    for (k = 0; k < n; k++) a[k] = x[k];
    if (typeof y !== "undefined")
        for (k = 0; k < n; k++) b[k] = y[k];
    cx[0] = 1;
    for (k = 1; k <= m / 2; k++) {
        t = c * k * k;
        cx[k] = cos(t);
        cy[k] = sin(t);
        cx[m - k] = cos(t);
        cy[m - k] = sin(t)
    }
    var X = new numeric.T(a, b),
        Y = new numeric.T(cx, cy);
    X = X.mul(Y);
    numeric.convpow2(X.x, X.y, numeric.clone(Y.x), numeric.neg(Y.y));
    X = X.mul(Y);
    X.x.length = n;
    X.y.length = n;
    return X
};
numeric.T.prototype.ifft = function ifft() {
    var x = this.x,
        y = this.y;
    var n = x.length,
        log = Math.log,
        log2 = log(2),
        p = Math.ceil(log(2 * n - 1) / log2),
        m = Math.pow(2, p);
    var cx = numeric.rep([m], 0),
        cy = numeric.rep([m], 0),
        cos = Math.cos,
        sin = Math.sin;
    var k, c = 3.141592653589793 / n,
        t;
    var a = numeric.rep([m], 0),
        b = numeric.rep([m], 0),
        nhalf = Math.floor(n / 2);
    for (k = 0; k < n; k++) a[k] = x[k];
    if (typeof y !== "undefined")
        for (k = 0; k < n; k++) b[k] = y[k];
    cx[0] = 1;
    for (k = 1; k <= m / 2; k++) {
        t = c * k * k;
        cx[k] = cos(t);
        cy[k] = sin(t);
        cx[m - k] = cos(t);
        cy[m - k] = sin(t)
    }
    var X = new numeric.T(a, b),
        Y = new numeric.T(cx, cy);
    X = X.mul(Y);
    numeric.convpow2(X.x, X.y, numeric.clone(Y.x), numeric.neg(Y.y));
    X = X.mul(Y);
    X.x.length = n;
    X.y.length = n;
    return X.div(n)
};
numeric.gradient = function gradient(f, x) {
    var n = x.length;
    var f0 = f(x);
    if (isNaN(f0)) throw new Error("gradient: f(x) is a NaN!");
    var max = Math.max;
    var i, x0 = numeric.clone(x),
        f1, f2, J = Array(n);
    var div = numeric.div,
        sub = numeric.sub,
        errest, roundoff, max = Math.max,
        eps = .001,
        abs = Math.abs,
        min = Math.min;
    var t0, t1, t2, it = 0,
        d1, d2, N;
    for (i = 0; i < n; i++) {
        var h = max(1e-6 * f0, 1e-8);
        while (1) {
            ++it;
            if (it > 20) {
                throw new Error("Numerical gradient fails")
            }
            x0[i] = x[i] + h;
            f1 = f(x0);
            x0[i] = x[i] - h;
            f2 = f(x0);
            x0[i] = x[i];
            if (isNaN(f1) || isNaN(f2)) {
                h /= 16;
                continue
            }
            J[i] = (f1 - f2) / (2 * h);
            t0 = x[i] - h;
            t1 = x[i];
            t2 = x[i] + h;
            d1 = (f1 - f0) / h;
            d2 = (f0 - f2) / h;
            N = max(abs(J[i]), abs(f0), abs(f1), abs(f2), abs(t0), abs(t1), abs(t2), 1e-8);
            errest = min(max(abs(d1 - J[i]), abs(d2 - J[i]), abs(d1 - d2)) / N, h / N);
            if (errest > eps) {
                h /= 16
            } else break
        }
    }
    return J
};
numeric.uncmin = function uncmin(f, x0, tol, gradient, maxit, callback, options) {
    var grad = numeric.gradient;
    if (typeof options === "undefined") {
        options = {}
    }
    if (typeof tol === "undefined") {
        tol = 1e-8
    }
    if (typeof gradient === "undefined") {
        gradient = function(x) {
            return grad(f, x)
        }
    }
    if (typeof maxit === "undefined") maxit = 1e3;
    x0 = numeric.clone(x0);
    var n = x0.length;
    var f0 = f(x0),
        f1, df0;
    if (isNaN(f0)) throw new Error("uncmin: f(x0) is a NaN!");
    var max = Math.max,
        norm2 = numeric.norm2;
    tol = max(tol, numeric.epsilon);
    var step, g0, g1, H1 = options.Hinv || numeric.identity(n);
    var dot = numeric.dot,
        inv = numeric.inv,
        sub = numeric.sub,
        add = numeric.add,
        ten = numeric.tensor,
        div = numeric.div,
        mul = numeric.mul;
    var all = numeric.all,
        isfinite = numeric.isFinite,
        neg = numeric.neg;
    var it = 0,
        i, s, x1, y, Hy, Hs, ys, i0, t, nstep, t1, t2;
    var msg = "";
    g0 = gradient(x0);
    while (it < maxit) {
        if (typeof callback === "function") {
            if (callback(it, x0, f0, g0, H1)) {
                msg = "Callback returned true";
                break
            }
        }
        if (!all(isfinite(g0))) {
            msg = "Gradient has Infinity or NaN";
            break
        }
        step = neg(dot(H1, g0));
        if (!all(isfinite(step))) {
            msg = "Search direction has Infinity or NaN";
            break
        }
        nstep = norm2(step);
        if (nstep < tol) {
            msg = "Newton step smaller than tol";
            break
        }
        t = 1;
        df0 = dot(g0, step);
        x1 = x0;
        while (it < maxit) {
            if (t * nstep < tol) {
                break
            }
            s = mul(step, t);
            x1 = add(x0, s);
            f1 = f(x1);
            if (f1 - f0 >= .1 * t * df0 || isNaN(f1)) {
                t *= .5;
                ++it;
                continue
            }
            break
        }
        if (t * nstep < tol) {
            msg = "Line search step size smaller than tol";
            break
        }
        if (it === maxit) {
            msg = "maxit reached during line search";
            break
        }
        g1 = gradient(x1);
        y = sub(g1, g0);
        ys = dot(y, s);
        Hy = dot(H1, y);
        H1 = sub(add(H1, mul((ys + dot(y, Hy)) / (ys * ys), ten(s, s))), div(add(ten(Hy, s), ten(s, Hy)), ys));
        x0 = x1;
        f0 = f1;
        g0 = g1;
        ++it
    }
    return {
        solution: x0,
        f: f0,
        gradient: g0,
        invHessian: H1,
        iterations: it,
        message: msg
    }
};
numeric.Dopri = function Dopri(x, y, f, ymid, iterations, msg, events) {
    this.x = x;
    this.y = y;
    this.f = f;
    this.ymid = ymid;
    this.iterations = iterations;
    this.events = events;
    this.message = msg
};
numeric.Dopri.prototype._at = function _at(xi, j) {
    function sqr(x) {
        return x * x
    }
    var sol = this;
    var xs = sol.x;
    var ys = sol.y;
    var k1 = sol.f;
    var ymid = sol.ymid;
    var n = xs.length;
    var x0, x1, xh, y0, y1, yh, xi;
    var floor = Math.floor,
        h;
    var c = .5;
    var add = numeric.add,
        mul = numeric.mul,
        sub = numeric.sub,
        p, q, w;
    x0 = xs[j];
    x1 = xs[j + 1];
    y0 = ys[j];
    y1 = ys[j + 1];
    h = x1 - x0;
    xh = x0 + c * h;
    yh = ymid[j];
    p = sub(k1[j], mul(y0, 1 / (x0 - xh) + 2 / (x0 - x1)));
    q = sub(k1[j + 1], mul(y1, 1 / (x1 - xh) + 2 / (x1 - x0)));
    w = [sqr(xi - x1) * (xi - xh) / sqr(x0 - x1) / (x0 - xh), sqr(xi - x0) * sqr(xi - x1) / sqr(x0 - xh) / sqr(x1 - xh), sqr(xi - x0) * (xi - xh) / sqr(x1 - x0) / (x1 - xh), (xi - x0) * sqr(xi - x1) * (xi - xh) / sqr(x0 - x1) / (x0 - xh), (xi - x1) * sqr(xi - x0) * (xi - xh) / sqr(x0 - x1) / (x1 - xh)];
    return add(add(add(add(mul(y0, w[0]), mul(yh, w[1])), mul(y1, w[2])), mul(p, w[3])), mul(q, w[4]))
};
numeric.Dopri.prototype.at = function at(x) {
    var i, j, k, floor = Math.floor;
    if (typeof x !== "number") {
        var n = x.length,
            ret = Array(n);
        for (i = n - 1; i !== -1; --i) {
            ret[i] = this.at(x[i])
        }
        return ret
    }
    var x0 = this.x;
    i = 0;
    j = x0.length - 1;
    while (j - i > 1) {
        k = floor(.5 * (i + j));
        if (x0[k] <= x) i = k;
        else j = k
    }
    return this._at(x, i)
};
numeric.dopri = function dopri(x0, x1, y0, f, tol, maxit, event) {
    if (typeof tol === "undefined") {
        tol = 1e-6
    }
    if (typeof maxit === "undefined") {
        maxit = 1e3
    }
    var xs = [x0],
        ys = [y0],
        k1 = [f(x0, y0)],
        k2, k3, k4, k5, k6, k7, ymid = [];
    var A2 = 1 / 5;
    var A3 = [3 / 40, 9 / 40];
    var A4 = [44 / 45, -56 / 15, 32 / 9];
    var A5 = [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729];
    var A6 = [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656];
    var b = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84];
    var bm = [.5 * 6025192743 / 30085553152, 0, .5 * 51252292925 / 65400821598, .5 * -2691868925 / 45128329728, .5 * 187940372067 / 1594534317056, .5 * -1776094331 / 19743644256, .5 * 11237099 / 235043384];
    var c = [1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
    var e = [-71 / 57600, 0, 71 / 16695, -71 / 1920, 17253 / 339200, -22 / 525, 1 / 40];
    var i = 0,
        er, j;
    var h = (x1 - x0) / 10;
    var it = 0;
    var add = numeric.add,
        mul = numeric.mul,
        y1, erinf;
    var max = Math.max,
        min = Math.min,
        abs = Math.abs,
        norminf = numeric.norminf,
        pow = Math.pow;
    var any = numeric.any,
        lt = numeric.lt,
        and = numeric.and,
        sub = numeric.sub;
    var e0, e1, ev;
    var ret = new numeric.Dopri(xs, ys, k1, ymid, -1, "");
    if (typeof event === "function") e0 = event(x0, y0);
    while (x0 < x1 && it < maxit) {
        ++it;
        if (x0 + h > x1) h = x1 - x0;
        k2 = f(x0 + c[0] * h, add(y0, mul(A2 * h, k1[i])));
        k3 = f(x0 + c[1] * h, add(add(y0, mul(A3[0] * h, k1[i])), mul(A3[1] * h, k2)));
        k4 = f(x0 + c[2] * h, add(add(add(y0, mul(A4[0] * h, k1[i])), mul(A4[1] * h, k2)), mul(A4[2] * h, k3)));
        k5 = f(x0 + c[3] * h, add(add(add(add(y0, mul(A5[0] * h, k1[i])), mul(A5[1] * h, k2)), mul(A5[2] * h, k3)), mul(A5[3] * h, k4)));
        k6 = f(x0 + c[4] * h, add(add(add(add(add(y0, mul(A6[0] * h, k1[i])), mul(A6[1] * h, k2)), mul(A6[2] * h, k3)), mul(A6[3] * h, k4)), mul(A6[4] * h, k5)));
        y1 = add(add(add(add(add(y0, mul(k1[i], h * b[0])), mul(k3, h * b[2])), mul(k4, h * b[3])), mul(k5, h * b[4])), mul(k6, h * b[5]));
        k7 = f(x0 + h, y1);
        er = add(add(add(add(add(mul(k1[i], h * e[0]), mul(k3, h * e[2])), mul(k4, h * e[3])), mul(k5, h * e[4])), mul(k6, h * e[5])), mul(k7, h * e[6]));
        if (typeof er === "number") erinf = abs(er);
        else erinf = norminf(er);
        if (erinf > tol) {
            h = .2 * h * pow(tol / erinf, .25);
            if (x0 + h === x0) {
                ret.msg = "Step size became too small";
                break
            }
            continue
        }
        ymid[i] = add(add(add(add(add(add(y0, mul(k1[i], h * bm[0])), mul(k3, h * bm[2])), mul(k4, h * bm[3])), mul(k5, h * bm[4])), mul(k6, h * bm[5])), mul(k7, h * bm[6]));
        ++i;
        xs[i] = x0 + h;
        ys[i] = y1;
        k1[i] = k7;
        if (typeof event === "function") {
            var yi, xl = x0,
                xr = x0 + .5 * h,
                xi;
            e1 = event(xr, ymid[i - 1]);
            ev = and(lt(e0, 0), lt(0, e1));
            if (!any(ev)) {
                xl = xr;
                xr = x0 + h;
                e0 = e1;
                e1 = event(xr, y1);
                ev = and(lt(e0, 0), lt(0, e1))
            }
            if (any(ev)) {
                var xc, yc, en, ei;
                var side = 0,
                    sl = 1,
                    sr = 1;
                while (1) {
                    if (typeof e0 === "number") xi = (sr * e1 * xl - sl * e0 * xr) / (sr * e1 - sl * e0);
                    else {
                        xi = xr;
                        for (j = e0.length - 1; j !== -1; --j) {
                            if (e0[j] < 0 && e1[j] > 0) xi = min(xi, (sr * e1[j] * xl - sl * e0[j] * xr) / (sr * e1[j] - sl * e0[j]))
                        }
                    }
                    if (xi <= xl || xi >= xr) break;
                    yi = ret._at(xi, i - 1);
                    ei = event(xi, yi);
                    en = and(lt(e0, 0), lt(0, ei));
                    if (any(en)) {
                        xr = xi;
                        e1 = ei;
                        ev = en;
                        sr = 1;
                        if (side === -1) sl *= .5;
                        else sl = 1;
                        side = -1
                    } else {
                        xl = xi;
                        e0 = ei;
                        sl = 1;
                        if (side === 1) sr *= .5;
                        else sr = 1;
                        side = 1
                    }
                }
                y1 = ret._at(.5 * (x0 + xi), i - 1);
                ret.f[i] = f(xi, yi);
                ret.x[i] = xi;
                ret.y[i] = yi;
                ret.ymid[i - 1] = y1;
                ret.events = ev;
                ret.iterations = it;
                return ret
            }
        }
        x0 += h;
        y0 = y1;
        e0 = e1;
        h = min(.8 * h * pow(tol / erinf, .25), 4 * h)
    }
    ret.iterations = it;
    return ret
};
numeric.LU = function(A, fast) {
    fast = fast || false;
    var abs = Math.abs;
    var i, j, k, absAjk, Akk, Ak, Pk, Ai;
    var max;
    var n = A.length,
        n1 = n - 1;
    var P = new Array(n);
    if (!fast) A = numeric.clone(A);
    for (k = 0; k < n; ++k) {
        Pk = k;
        Ak = A[k];
        max = abs(Ak[k]);
        for (j = k + 1; j < n; ++j) {
            absAjk = abs(A[j][k]);
            if (max < absAjk) {
                max = absAjk;
                Pk = j
            }
        }
        P[k] = Pk;
        if (Pk != k) {
            A[k] = A[Pk];
            A[Pk] = Ak;
            Ak = A[k]
        }
        Akk = Ak[k];
        for (i = k + 1; i < n; ++i) {
            A[i][k] /= Akk
        }
        for (i = k + 1; i < n; ++i) {
            Ai = A[i];
            for (j = k + 1; j < n1; ++j) {
                Ai[j] -= Ai[k] * Ak[j];
                ++j;
                Ai[j] -= Ai[k] * Ak[j]
            }
            if (j === n1) Ai[j] -= Ai[k] * Ak[j]
        }
    }
    return {
        LU: A,
        P: P
    }
};
numeric.LUsolve = function LUsolve(LUP, b) {
    var i, j;
    var LU = LUP.LU;
    var n = LU.length;
    var x = numeric.clone(b);
    var P = LUP.P;
    var Pi, LUi, LUii, tmp;
    for (i = n - 1; i !== -1; --i) x[i] = b[i];
    for (i = 0; i < n; ++i) {
        Pi = P[i];
        if (P[i] !== i) {
            tmp = x[i];
            x[i] = x[Pi];
            x[Pi] = tmp
        }
        LUi = LU[i];
        for (j = 0; j < i; ++j) {
            x[i] -= x[j] * LUi[j]
        }
    }
    for (i = n - 1; i >= 0; --i) {
        LUi = LU[i];
        for (j = i + 1; j < n; ++j) {
            x[i] -= x[j] * LUi[j]
        }
        x[i] /= LUi[i]
    }
    return x
};
numeric.solve = function solve(A, b, fast) {
    return numeric.LUsolve(numeric.LU(A, fast), b)
};
numeric.echelonize = function echelonize(A) {
    var s = numeric.dim(A),
        m = s[0],
        n = s[1];
    var I = numeric.identity(m);
    var P = Array(m);
    var i, j, k, l, Ai, Ii, Z, a;
    var abs = Math.abs;
    var diveq = numeric.diveq;
    A = numeric.clone(A);
    for (i = 0; i < m; ++i) {
        k = 0;
        Ai = A[i];
        Ii = I[i];
        for (j = 1; j < n; ++j)
            if (abs(Ai[k]) < abs(Ai[j])) k = j;
        P[i] = k;
        diveq(Ii, Ai[k]);
        diveq(Ai, Ai[k]);
        for (j = 0; j < m; ++j)
            if (j !== i) {
                Z = A[j];
                a = Z[k];
                for (l = n - 1; l !== -1; --l) Z[l] -= Ai[l] * a;
                Z = I[j];
                for (l = m - 1; l !== -1; --l) Z[l] -= Ii[l] * a
            }
    }
    return {
        I: I,
        A: A,
        P: P
    }
};
numeric.__solveLP = function __solveLP(c, A, b, tol, maxit, x, flag) {
    var sum = numeric.sum,
        log = numeric.log,
        mul = numeric.mul,
        sub = numeric.sub,
        dot = numeric.dot,
        div = numeric.div,
        add = numeric.add;
    var m = c.length,
        n = b.length,
        y;
    var unbounded = false,
        cb, i0 = 0;
    var alpha = 1;
    var f0, df0, AT = numeric.transpose(A),
        svd = numeric.svd,
        transpose = numeric.transpose,
        leq = numeric.leq,
        sqrt = Math.sqrt,
        abs = Math.abs;
    var muleq = numeric.muleq;
    var norm = numeric.norminf,
        any = numeric.any,
        min = Math.min;
    var all = numeric.all,
        gt = numeric.gt;
    var p = Array(m),
        A0 = Array(n),
        e = numeric.rep([n], 1),
        H;
    var solve = numeric.solve,
        z = sub(b, dot(A, x)),
        count;
    var dotcc = dot(c, c);
    var g;
    for (count = i0; count < maxit; ++count) {
        var i, j, d;
        for (i = n - 1; i !== -1; --i) A0[i] = div(A[i], z[i]);
        var A1 = transpose(A0);
        for (i = m - 1; i !== -1; --i) p[i] = sum(A1[i]);
        alpha = .25 * abs(dotcc / dot(c, p));
        var a1 = 100 * sqrt(dotcc / dot(p, p));
        if (!isFinite(alpha) || alpha > a1) alpha = a1;
        g = add(c, mul(alpha, p));
        H = dot(A1, A0);
        for (i = m - 1; i !== -1; --i) H[i][i] += 1;
        d = solve(H, div(g, alpha), true);
        var t0 = div(z, dot(A, d));
        var t = 1;
        for (i = n - 1; i !== -1; --i)
            if (t0[i] < 0) t = min(t, -.999 * t0[i]);
        y = sub(x, mul(d, t));
        z = sub(b, dot(A, y));
        if (!all(gt(z, 0))) return {
            solution: x,
            message: "",
            iterations: count
        };
        x = y;
        if (alpha < tol) return {
            solution: y,
            message: "",
            iterations: count
        };
        if (flag) {
            var s = dot(c, g),
                Ag = dot(A, g);
            unbounded = true;
            for (i = n - 1; i !== -1; --i)
                if (s * Ag[i] < 0) {
                    unbounded = false;
                    break
                }
        } else {
            if (x[m - 1] >= 0) unbounded = false;
            else unbounded = true
        }
        if (unbounded) return {
            solution: y,
            message: "Unbounded",
            iterations: count
        }
    }
    return {
        solution: x,
        message: "maximum iteration count exceeded",
        iterations: count
    }
};
numeric._solveLP = function _solveLP(c, A, b, tol, maxit) {
    var m = c.length,
        n = b.length,
        y;
    var sum = numeric.sum,
        log = numeric.log,
        mul = numeric.mul,
        sub = numeric.sub,
        dot = numeric.dot,
        div = numeric.div,
        add = numeric.add;
    var c0 = numeric.rep([m], 0).concat([1]);
    var J = numeric.rep([n, 1], -1);
    var A0 = numeric.blockMatrix([
        [A, J]
    ]);
    var b0 = b;
    var y = numeric.rep([m], 0).concat(Math.max(0, numeric.sup(numeric.neg(b))) + 1);
    var x0 = numeric.__solveLP(c0, A0, b0, tol, maxit, y, false);
    var x = numeric.clone(x0.solution);
    x.length = m;
    var foo = numeric.inf(sub(b, dot(A, x)));
    if (foo < 0) {
        return {
            solution: NaN,
            message: "Infeasible",
            iterations: x0.iterations
        }
    }
    var ret = numeric.__solveLP(c, A, b, tol, maxit - x0.iterations, x, true);
    ret.iterations += x0.iterations;
    return ret
};
numeric.solveLP = function solveLP(c, A, b, Aeq, beq, tol, maxit) {
    if (typeof maxit === "undefined") maxit = 1e3;
    if (typeof tol === "undefined") tol = numeric.epsilon;
    if (typeof Aeq === "undefined") return numeric._solveLP(c, A, b, tol, maxit);
    var m = Aeq.length,
        n = Aeq[0].length,
        o = A.length;
    var B = numeric.echelonize(Aeq);
    var flags = numeric.rep([n], 0);
    var P = B.P;
    var Q = [];
    var i;
    for (i = P.length - 1; i !== -1; --i) flags[P[i]] = 1;
    for (i = n - 1; i !== -1; --i)
        if (flags[i] === 0) Q.push(i);
    var g = numeric.getRange;
    var I = numeric.linspace(0, m - 1),
        J = numeric.linspace(0, o - 1);
    var Aeq2 = g(Aeq, I, Q),
        A1 = g(A, J, P),
        A2 = g(A, J, Q),
        dot = numeric.dot,
        sub = numeric.sub;
    var A3 = dot(A1, B.I);
    var A4 = sub(A2, dot(A3, Aeq2)),
        b4 = sub(b, dot(A3, beq));
    var c1 = Array(P.length),
        c2 = Array(Q.length);
    for (i = P.length - 1; i !== -1; --i) c1[i] = c[P[i]];
    for (i = Q.length - 1; i !== -1; --i) c2[i] = c[Q[i]];
    var c4 = sub(c2, dot(c1, dot(B.I, Aeq2)));
    var S = numeric._solveLP(c4, A4, b4, tol, maxit);
    var x2 = S.solution;
    if (x2 !== x2) return S;
    var x1 = dot(B.I, sub(beq, dot(Aeq2, x2)));
    var x = Array(c.length);
    for (i = P.length - 1; i !== -1; --i) x[P[i]] = x1[i];
    for (i = Q.length - 1; i !== -1; --i) x[Q[i]] = x2[i];
    return {
        solution: x,
        message: S.message,
        iterations: S.iterations
    }
};
numeric.MPStoLP = function MPStoLP(MPS) {
    if (MPS instanceof String) {
        MPS.split("\n")
    }
    var state = 0;
    var states = ["Initial state", "NAME", "ROWS", "COLUMNS", "RHS", "BOUNDS", "ENDATA"];
    var n = MPS.length;
    var i, j, z, N = 0,
        rows = {},
        sign = [],
        rl = 0,
        vars = {},
        nv = 0;
    var name;
    var c = [],
        A = [],
        b = [];

    function err(e) {
        throw new Error("MPStoLP: " + e + "\nLine " + i + ": " + MPS[i] + "\nCurrent state: " + states[state] + "\n")
    }
    for (i = 0; i < n; ++i) {
        z = MPS[i];
        var w0 = z.match(/\S*/g);
        var w = [];
        for (j = 0; j < w0.length; ++j)
            if (w0[j] !== "") w.push(w0[j]);
        if (w.length === 0) continue;
        for (j = 0; j < states.length; ++j)
            if (z.substr(0, states[j].length) === states[j]) break;
        if (j < states.length) {
            state = j;
            if (j === 1) {
                name = w[1]
            }
            if (j === 6) return {
                name: name,
                c: c,
                A: numeric.transpose(A),
                b: b,
                rows: rows,
                vars: vars
            };
            continue
        }
        switch (state) {
            case 0:
            case 1:
                err("Unexpected line");
            case 2:
                switch (w[0]) {
                    case "N":
                        if (N === 0) N = w[1];
                        else err("Two or more N rows");
                        break;
                    case "L":
                        rows[w[1]] = rl;
                        sign[rl] = 1;
                        b[rl] = 0;
                        ++rl;
                        break;
                    case "G":
                        rows[w[1]] = rl;
                        sign[rl] = -1;
                        b[rl] = 0;
                        ++rl;
                        break;
                    case "E":
                        rows[w[1]] = rl;
                        sign[rl] = 0;
                        b[rl] = 0;
                        ++rl;
                        break;
                    default:
                        err("Parse error " + numeric.prettyPrint(w))
                }
                break;
            case 3:
                if (!vars.hasOwnProperty(w[0])) {
                    vars[w[0]] = nv;
                    c[nv] = 0;
                    A[nv] = numeric.rep([rl], 0);
                    ++nv
                }
                var p = vars[w[0]];
                for (j = 1; j < w.length; j += 2) {
                    if (w[j] === N) {
                        c[p] = parseFloat(w[j + 1]);
                        continue
                    }
                    var q = rows[w[j]];
                    A[p][q] = (sign[q] < 0 ? -1 : 1) * parseFloat(w[j + 1])
                }
                break;
            case 4:
                for (j = 1; j < w.length; j += 2) b[rows[w[j]]] = (sign[rows[w[j]]] < 0 ? -1 : 1) * parseFloat(w[j + 1]);
                break;
            case 5:
                break;
            case 6:
                err("Internal error")
        }
    }
    err("Reached end of file without ENDATA")
};
numeric.seedrandom = {
    pow: Math.pow,
    random: Math.random
};
(function(pool, math, width, chunks, significance, overflow, startdenom) {
    math["seedrandom"] = function seedrandom(seed, use_entropy) {
        var key = [];
        var arc4;
        seed = mixkey(flatten(use_entropy ? [seed, pool] : arguments.length ? seed : [(new Date).getTime(), pool, window], 3), key);
        arc4 = new ARC4(key);
        mixkey(arc4.S, pool);
        math["random"] = function random() {
            var n = arc4.g(chunks);
            var d = startdenom;
            var x = 0;
            while (n < significance) {
                n = (n + x) * width;
                d *= width;
                x = arc4.g(1)
            }
            while (n >= overflow) {
                n /= 2;
                d /= 2;
                x >>>= 1
            }
            return (n + x) / d
        };
        return seed
    };

    function ARC4(key) {
        var t, u, me = this,
            keylen = key.length;
        var i = 0,
            j = me.i = me.j = me.m = 0;
        me.S = [];
        me.c = [];
        if (!keylen) {
            key = [keylen++]
        }
        while (i < width) {
            me.S[i] = i++
        }
        for (i = 0; i < width; i++) {
            t = me.S[i];
            j = lowbits(j + t + key[i % keylen]);
            u = me.S[j];
            me.S[i] = u;
            me.S[j] = t
        }
        me.g = function getnext(count) {
            var s = me.S;
            var i = lowbits(me.i + 1);
            var t = s[i];
            var j = lowbits(me.j + t);
            var u = s[j];
            s[i] = u;
            s[j] = t;
            var r = s[lowbits(t + u)];
            while (--count) {
                i = lowbits(i + 1);
                t = s[i];
                j = lowbits(j + t);
                u = s[j];
                s[i] = u;
                s[j] = t;
                r = r * width + s[lowbits(t + u)]
            }
            me.i = i;
            me.j = j;
            return r
        };
        me.g(width)
    }

    function flatten(obj, depth, result, prop, typ) {
        result = [];
        typ = typeof obj;
        if (depth && typ == "object") {
            for (prop in obj) {
                if (prop.indexOf("S") < 5) {
                    try {
                        result.push(flatten(obj[prop], depth - 1))
                    } catch (e) {}
                }
            }
        }
        return result.length ? result : obj + (typ != "string" ? "\x00" : "")
    }

    function mixkey(seed, key, smear, j) {
        seed += "";
        smear = 0;
        for (j = 0; j < seed.length; j++) {
            key[lowbits(j)] = lowbits((smear ^= key[lowbits(j)] * 19) + seed.charCodeAt(j))
        }
        seed = "";
        for (j in key) {
            seed += String.fromCharCode(key[j])
        }
        return seed
    }

    function lowbits(n) {
        return n & width - 1
    }
    startdenom = math.pow(width, chunks);
    significance = math.pow(2, significance);
    overflow = significance * 2;
    mixkey(math.random(), pool)
})([], numeric.seedrandom, 256, 6, 52);
(function(exports) {
    function base0to1(A) {
        if (typeof A !== "object") {
            return A
        }
        var ret = [],
            i, n = A.length;
        for (i = 0; i < n; i++) ret[i + 1] = base0to1(A[i]);
        return ret
    }

    function base1to0(A) {
        if (typeof A !== "object") {
            return A
        }
        var ret = [],
            i, n = A.length;
        for (i = 1; i < n; i++) ret[i - 1] = base1to0(A[i]);
        return ret
    }

    function dpori(a, lda, n) {
        var i, j, k, kp1, t;
        for (k = 1; k <= n; k = k + 1) {
            a[k][k] = 1 / a[k][k];
            t = -a[k][k];
            for (i = 1; i < k; i = i + 1) {
                a[i][k] = t * a[i][k]
            }
            kp1 = k + 1;
            if (n < kp1) {
                break
            }
            for (j = kp1; j <= n; j = j + 1) {
                t = a[k][j];
                a[k][j] = 0;
                for (i = 1; i <= k; i = i + 1) {
                    a[i][j] = a[i][j] + t * a[i][k]
                }
            }
        }
    }

    function dposl(a, lda, n, b) {
        var i, k, kb, t;
        for (k = 1; k <= n; k = k + 1) {
            t = 0;
            for (i = 1; i < k; i = i + 1) {
                t = t + a[i][k] * b[i]
            }
            b[k] = (b[k] - t) / a[k][k]
        }
        for (kb = 1; kb <= n; kb = kb + 1) {
            k = n + 1 - kb;
            b[k] = b[k] / a[k][k];
            t = -b[k];
            for (i = 1; i < k; i = i + 1) {
                b[i] = b[i] + t * a[i][k]
            }
        }
    }

    function dpofa(a, lda, n, info) {
        var i, j, jm1, k, t, s;
        for (j = 1; j <= n; j = j + 1) {
            info[1] = j;
            s = 0;
            jm1 = j - 1;
            if (jm1 < 1) {
                s = a[j][j] - s;
                if (s <= 0) {
                    break
                }
                a[j][j] = Math.sqrt(s)
            } else {
                for (k = 1; k <= jm1; k = k + 1) {
                    t = a[k][j];
                    for (i = 1; i < k; i = i + 1) {
                        t = t - a[i][j] * a[i][k]
                    }
                    t = t / a[k][k];
                    a[k][j] = t;
                    s = s + t * t
                }
                s = a[j][j] - s;
                if (s <= 0) {
                    break
                }
                a[j][j] = Math.sqrt(s)
            }
            info[1] = 0
        }
    }

    function qpgen2(dmat, dvec, fddmat, n, sol, crval, amat, bvec, fdamat, q, meq, iact, nact, iter, work, ierr) {
        var i, j, l, l1, info, it1, iwzv, iwrv, iwrm, iwsv, iwuv, nvl, r, iwnbv, temp, sum, t1, tt, gc, gs, nu, t1inf, t2min, vsmall, tmpa, tmpb, go;
        r = Math.min(n, q);
        l = 2 * n + r * (r + 5) / 2 + 2 * q + 1;
        vsmall = 1e-60;
        do {
            vsmall = vsmall + vsmall;
            tmpa = 1 + .1 * vsmall;
            tmpb = 1 + .2 * vsmall
        } while (tmpa <= 1 || tmpb <= 1);
        for (i = 1; i <= n; i = i + 1) {
            work[i] = dvec[i]
        }
        for (i = n + 1; i <= l; i = i + 1) {
            work[i] = 0
        }
        for (i = 1; i <= q; i = i + 1) {
            iact[i] = 0
        }
        info = [];
        if (ierr[1] === 0) {
            dpofa(dmat, fddmat, n, info);
            if (info[1] !== 0) {
                ierr[1] = 2;
                return
            }
            dposl(dmat, fddmat, n, dvec);
            dpori(dmat, fddmat, n)
        } else {
            for (j = 1; j <= n; j = j + 1) {
                sol[j] = 0;
                for (i = 1; i <= j; i = i + 1) {
                    sol[j] = sol[j] + dmat[i][j] * dvec[i]
                }
            }
            for (j = 1; j <= n; j = j + 1) {
                dvec[j] = 0;
                for (i = j; i <= n; i = i + 1) {
                    dvec[j] = dvec[j] + dmat[j][i] * sol[i]
                }
            }
        }
        crval[1] = 0;
        for (j = 1; j <= n; j = j + 1) {
            sol[j] = dvec[j];
            crval[1] = crval[1] + work[j] * sol[j];
            work[j] = 0;
            for (i = j + 1; i <= n; i = i + 1) {
                dmat[i][j] = 0
            }
        }
        crval[1] = -crval[1] / 2;
        ierr[1] = 0;
        iwzv = n;
        iwrv = iwzv + n;
        iwuv = iwrv + r;
        iwrm = iwuv + r + 1;
        iwsv = iwrm + r * (r + 1) / 2;
        iwnbv = iwsv + q;
        for (i = 1; i <= q; i = i + 1) {
            sum = 0;
            for (j = 1; j <= n; j = j + 1) {
                sum = sum + amat[j][i] * amat[j][i]
            }
            work[iwnbv + i] = Math.sqrt(sum)
        }
        nact = 0;
        iter[1] = 0;
        iter[2] = 0;

        function fn_goto_50() {
            iter[1] = iter[1] + 1;
            l = iwsv;
            for (i = 1; i <= q; i = i + 1) {
                l = l + 1;
                sum = -bvec[i];
                for (j = 1; j <= n; j = j + 1) {
                    sum = sum + amat[j][i] * sol[j]
                }
                if (Math.abs(sum) < vsmall) {
                    sum = 0
                }
                if (i > meq) {
                    work[l] = sum
                } else {
                    work[l] = -Math.abs(sum);
                    if (sum > 0) {
                        for (j = 1; j <= n; j = j + 1) {
                            amat[j][i] = -amat[j][i]
                        }
                        bvec[i] = -bvec[i]
                    }
                }
            }
            for (i = 1; i <= nact; i = i + 1) {
                work[iwsv + iact[i]] = 0
            }
            nvl = 0;
            temp = 0;
            for (i = 1; i <= q; i = i + 1) {
                if (work[iwsv + i] < temp * work[iwnbv + i]) {
                    nvl = i;
                    temp = work[iwsv + i] / work[iwnbv + i]
                }
            }
            if (nvl === 0) {
                return 999
            }
            return 0
        }

        function fn_goto_55() {
            for (i = 1; i <= n; i = i + 1) {
                sum = 0;
                for (j = 1; j <= n; j = j + 1) {
                    sum = sum + dmat[j][i] * amat[j][nvl]
                }
                work[i] = sum
            }
            l1 = iwzv;
            for (i = 1; i <= n; i = i + 1) {
                work[l1 + i] = 0
            }
            for (j = nact + 1; j <= n; j = j + 1) {
                for (i = 1; i <= n; i = i + 1) {
                    work[l1 + i] = work[l1 + i] + dmat[i][j] * work[j]
                }
            }
            t1inf = true;
            for (i = nact; i >= 1; i = i - 1) {
                sum = work[i];
                l = iwrm + i * (i + 3) / 2;
                l1 = l - i;
                for (j = i + 1; j <= nact; j = j + 1) {
                    sum = sum - work[l] * work[iwrv + j];
                    l = l + j
                }
                sum = sum / work[l1];
                work[iwrv + i] = sum;
                if (iact[i] < meq) {
                    break
                }
                if (sum < 0) {
                    break
                }
                t1inf = false;
                it1 = i
            }
            if (!t1inf) {
                t1 = work[iwuv + it1] / work[iwrv + it1];
                for (i = 1; i <= nact; i = i + 1) {
                    if (iact[i] < meq) {
                        break
                    }
                    if (work[iwrv + i] < 0) {
                        break
                    }
                    temp = work[iwuv + i] / work[iwrv + i];
                    if (temp < t1) {
                        t1 = temp;
                        it1 = i
                    }
                }
            }
            sum = 0;
            for (i = iwzv + 1; i <= iwzv + n; i = i + 1) {
                sum = sum + work[i] * work[i]
            }
            if (Math.abs(sum) <= vsmall) {
                if (t1inf) {
                    ierr[1] = 1;
                    return 999
                } else {
                    for (i = 1; i <= nact; i = i + 1) {
                        work[iwuv + i] = work[iwuv + i] - t1 * work[iwrv + i]
                    }
                    work[iwuv + nact + 1] = work[iwuv + nact + 1] + t1;
                    return 700
                }
            } else {
                sum = 0;
                for (i = 1; i <= n; i = i + 1) {
                    sum = sum + work[iwzv + i] * amat[i][nvl]
                }
                tt = -work[iwsv + nvl] / sum;
                t2min = true;
                if (!t1inf) {
                    if (t1 < tt) {
                        tt = t1;
                        t2min = false
                    }
                }
                for (i = 1; i <= n; i = i + 1) {
                    sol[i] = sol[i] + tt * work[iwzv + i];
                    if (Math.abs(sol[i]) < vsmall) {
                        sol[i] = 0
                    }
                }
                crval[1] = crval[1] + tt * sum * (tt / 2 + work[iwuv + nact + 1]);
                for (i = 1; i <= nact; i = i + 1) {
                    work[iwuv + i] = work[iwuv + i] - tt * work[iwrv + i]
                }
                work[iwuv + nact + 1] = work[iwuv + nact + 1] + tt;
                if (t2min) {
                    nact = nact + 1;
                    iact[nact] = nvl;
                    l = iwrm + (nact - 1) * nact / 2 + 1;
                    for (i = 1; i <= nact - 1; i = i + 1) {
                        work[l] = work[i];
                        l = l + 1
                    }
                    if (nact === n) {
                        work[l] = work[n]
                    } else {
                        for (i = n; i >= nact + 1; i = i - 1) {
                            if (work[i] === 0) {
                                break
                            }
                            gc = Math.max(Math.abs(work[i - 1]), Math.abs(work[i]));
                            gs = Math.min(Math.abs(work[i - 1]), Math.abs(work[i]));
                            if (work[i - 1] >= 0) {
                                temp = Math.abs(gc * Math.sqrt(1 + gs * gs / (gc * gc)))
                            } else {
                                temp = -Math.abs(gc * Math.sqrt(1 + gs * gs / (gc * gc)))
                            }
                            gc = work[i - 1] / temp;
                            gs = work[i] / temp;
                            if (gc === 1) {
                                break
                            }
                            if (gc === 0) {
                                work[i - 1] = gs * temp;
                                for (j = 1; j <= n; j = j + 1) {
                                    temp = dmat[j][i - 1];
                                    dmat[j][i - 1] = dmat[j][i];
                                    dmat[j][i] = temp
                                }
                            } else {
                                work[i - 1] = temp;
                                nu = gs / (1 + gc);
                                for (j = 1; j <= n; j = j + 1) {
                                    temp = gc * dmat[j][i - 1] + gs * dmat[j][i];
                                    dmat[j][i] = nu * (dmat[j][i - 1] + temp) - dmat[j][i];
                                    dmat[j][i - 1] = temp
                                }
                            }
                        }
                        work[l] = work[nact]
                    }
                } else {
                    sum = -bvec[nvl];
                    for (j = 1; j <= n; j = j + 1) {
                        sum = sum + sol[j] * amat[j][nvl]
                    }
                    if (nvl > meq) {
                        work[iwsv + nvl] = sum
                    } else {
                        work[iwsv + nvl] = -Math.abs(sum);
                        if (sum > 0) {
                            for (j = 1; j <= n; j = j + 1) {
                                amat[j][nvl] = -amat[j][nvl]
                            }
                            bvec[nvl] = -bvec[nvl]
                        }
                    }
                    return 700
                }
            }
            return 0
        }

        function fn_goto_797() {
            l = iwrm + it1 * (it1 + 1) / 2 + 1;
            l1 = l + it1;
            if (work[l1] === 0) {
                return 798
            }
            gc = Math.max(Math.abs(work[l1 - 1]), Math.abs(work[l1]));
            gs = Math.min(Math.abs(work[l1 - 1]), Math.abs(work[l1]));
            if (work[l1 - 1] >= 0) {
                temp = Math.abs(gc * Math.sqrt(1 + gs * gs / (gc * gc)))
            } else {
                temp = -Math.abs(gc * Math.sqrt(1 + gs * gs / (gc * gc)))
            }
            gc = work[l1 - 1] / temp;
            gs = work[l1] / temp;
            if (gc === 1) {
                return 798
            }
            if (gc === 0) {
                for (i = it1 + 1; i <= nact; i = i + 1) {
                    temp = work[l1 - 1];
                    work[l1 - 1] = work[l1];
                    work[l1] = temp;
                    l1 = l1 + i
                }
                for (i = 1; i <= n; i = i + 1) {
                    temp = dmat[i][it1];
                    dmat[i][it1] = dmat[i][it1 + 1];
                    dmat[i][it1 + 1] = temp
                }
            } else {
                nu = gs / (1 + gc);
                for (i = it1 + 1; i <= nact; i = i + 1) {
                    temp = gc * work[l1 - 1] + gs * work[l1];
                    work[l1] = nu * (work[l1 - 1] + temp) - work[l1];
                    work[l1 - 1] = temp;
                    l1 = l1 + i
                }
                for (i = 1; i <= n; i = i + 1) {
                    temp = gc * dmat[i][it1] + gs * dmat[i][it1 + 1];
                    dmat[i][it1 + 1] = nu * (dmat[i][it1] + temp) - dmat[i][it1 + 1];
                    dmat[i][it1] = temp
                }
            }
            return 0
        }

        function fn_goto_798() {
            l1 = l - it1;
            for (i = 1; i <= it1; i = i + 1) {
                work[l1] = work[l];
                l = l + 1;
                l1 = l1 + 1
            }
            work[iwuv + it1] = work[iwuv + it1 + 1];
            iact[it1] = iact[it1 + 1];
            it1 = it1 + 1;
            if (it1 < nact) {
                return 797
            }
            return 0
        }

        function fn_goto_799() {
            work[iwuv + nact] = work[iwuv + nact + 1];
            work[iwuv + nact + 1] = 0;
            iact[nact] = 0;
            nact = nact - 1;
            iter[2] = iter[2] + 1;
            return 0
        }
        go = 0;
        while (true) {
            go = fn_goto_50();
            if (go === 999) {
                return
            }
            while (true) {
                go = fn_goto_55();
                if (go === 0) {
                    break
                }
                if (go === 999) {
                    return
                }
                if (go === 700) {
                    if (it1 === nact) {
                        fn_goto_799()
                    } else {
                        while (true) {
                            fn_goto_797();
                            go = fn_goto_798();
                            if (go !== 797) {
                                break
                            }
                        }
                        fn_goto_799()
                    }
                }
            }
        }
    }

    function solveQP(Dmat, dvec, Amat, bvec, meq, factorized) {
        Dmat = base0to1(Dmat);
        dvec = base0to1(dvec);
        Amat = base0to1(Amat);
        var i, n, q, nact, r, crval = [],
            iact = [],
            sol = [],
            work = [],
            iter = [],
            message;
        meq = meq || 0;
        factorized = factorized ? base0to1(factorized) : [undefined, 0];
        bvec = bvec ? base0to1(bvec) : [];
        n = Dmat.length - 1;
        q = Amat[1].length - 1;
        if (!bvec) {
            for (i = 1; i <= q; i = i + 1) {
                bvec[i] = 0
            }
        }
        for (i = 1; i <= q; i = i + 1) {
            iact[i] = 0
        }
        nact = 0;
        r = Math.min(n, q);
        for (i = 1; i <= n; i = i + 1) {
            sol[i] = 0
        }
        crval[1] = 0;
        for (i = 1; i <= 2 * n + r * (r + 5) / 2 + 2 * q + 1; i = i + 1) {
            work[i] = 0
        }
        for (i = 1; i <= 2; i = i + 1) {
            iter[i] = 0
        }
        qpgen2(Dmat, dvec, n, n, sol, crval, Amat, bvec, n, q, meq, iact, nact, iter, work, factorized);
        message = "";
        if (factorized[1] === 1) {
            message = "constraints are inconsistent, no solution!"
        }
        if (factorized[1] === 2) {
            message = "matrix D in quadratic function is not positive definite!"
        }
        return {
            solution: base1to0(sol),
            value: base1to0(crval),
            unconstrained_solution: base1to0(dvec),
            iterations: base1to0(iter),
            iact: base1to0(iact),
            message: message
        }
    }
    exports.solveQP = solveQP
})(numeric);
numeric.svd = function svd(A) {
    var temp;
    var prec = numeric.epsilon;
    var tolerance = 1e-64 / prec;
    var itmax = 50;
    var c = 0;
    var i = 0;
    var j = 0;
    var k = 0;
    var l = 0;
    var u = numeric.clone(A);
    var m = u.length;
    var n = u[0].length;
    if (m < n) throw "Need more rows than columns";
    var e = new Array(n);
    var q = new Array(n);
    for (i = 0; i < n; i++) e[i] = q[i] = 0;
    var v = numeric.rep([n, n], 0);

    function pythag(a, b) {
        a = Math.abs(a);
        b = Math.abs(b);
        if (a > b) return a * Math.sqrt(1 + b * b / a / a);
        else if (b == 0) return a;
        return b * Math.sqrt(1 + a * a / b / b)
    }
    var f = 0;
    var g = 0;
    var h = 0;
    var x = 0;
    var y = 0;
    var z = 0;
    var s = 0;
    for (i = 0; i < n; i++) {
        e[i] = g;
        s = 0;
        l = i + 1;
        for (j = i; j < m; j++) s += u[j][i] * u[j][i];
        if (s <= tolerance) g = 0;
        else {
            f = u[i][i];
            g = Math.sqrt(s);
            if (f >= 0) g = -g;
            h = f * g - s;
            u[i][i] = f - g;
            for (j = l; j < n; j++) {
                s = 0;
                for (k = i; k < m; k++) s += u[k][i] * u[k][j];
                f = s / h;
                for (k = i; k < m; k++) u[k][j] += f * u[k][i]
            }
        }
        q[i] = g;
        s = 0;
        for (j = l; j < n; j++) s = s + u[i][j] * u[i][j];
        if (s <= tolerance) g = 0;
        else {
            f = u[i][i + 1];
            g = Math.sqrt(s);
            if (f >= 0) g = -g;
            h = f * g - s;
            u[i][i + 1] = f - g;
            for (j = l; j < n; j++) e[j] = u[i][j] / h;
            for (j = l; j < m; j++) {
                s = 0;
                for (k = l; k < n; k++) s += u[j][k] * u[i][k];
                for (k = l; k < n; k++) u[j][k] += s * e[k]
            }
        }
        y = Math.abs(q[i]) + Math.abs(e[i]);
        if (y > x) x = y
    }
    for (i = n - 1; i != -1; i += -1) {
        if (g != 0) {
            h = g * u[i][i + 1];
            for (j = l; j < n; j++) v[j][i] = u[i][j] / h;
            for (j = l; j < n; j++) {
                s = 0;
                for (k = l; k < n; k++) s += u[i][k] * v[k][j];
                for (k = l; k < n; k++) v[k][j] += s * v[k][i]
            }
        }
        for (j = l; j < n; j++) {
            v[i][j] = 0;
            v[j][i] = 0
        }
        v[i][i] = 1;
        g = e[i];
        l = i
    }
    for (i = n - 1; i != -1; i += -1) {
        l = i + 1;
        g = q[i];
        for (j = l; j < n; j++) u[i][j] = 0;
        if (g != 0) {
            h = u[i][i] * g;
            for (j = l; j < n; j++) {
                s = 0;
                for (k = l; k < m; k++) s += u[k][i] * u[k][j];
                f = s / h;
                for (k = i; k < m; k++) u[k][j] += f * u[k][i]
            }
            for (j = i; j < m; j++) u[j][i] = u[j][i] / g
        } else
            for (j = i; j < m; j++) u[j][i] = 0;
        u[i][i] += 1
    }
    prec = prec * x;
    for (k = n - 1; k != -1; k += -1) {
        for (var iteration = 0; iteration < itmax; iteration++) {
            var test_convergence = false;
            for (l = k; l != -1; l += -1) {
                if (Math.abs(e[l]) <= prec) {
                    test_convergence = true;
                    break
                }
                if (Math.abs(q[l - 1]) <= prec) break
            }
            if (!test_convergence) {
                c = 0;
                s = 1;
                var l1 = l - 1;
                for (i = l; i < k + 1; i++) {
                    f = s * e[i];
                    e[i] = c * e[i];
                    if (Math.abs(f) <= prec) break;
                    g = q[i];
                    h = pythag(f, g);
                    q[i] = h;
                    c = g / h;
                    s = -f / h;
                    for (j = 0; j < m; j++) {
                        y = u[j][l1];
                        z = u[j][i];
                        u[j][l1] = y * c + z * s;
                        u[j][i] = -y * s + z * c
                    }
                }
            }
            z = q[k];
            if (l == k) {
                if (z < 0) {
                    q[k] = -z;
                    for (j = 0; j < n; j++) v[j][k] = -v[j][k]
                }
                break
            }
            if (iteration >= itmax - 1) throw "Error: no convergence.";
            x = q[l];
            y = q[k - 1];
            g = e[k - 1];
            h = e[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
            g = pythag(f, 1);
            if (f < 0) f = ((x - z) * (x + z) + h * (y / (f - g) - h)) / x;
            else f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
            c = 1;
            s = 1;
            for (i = l + 1; i < k + 1; i++) {
                g = e[i];
                y = q[i];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                e[i - 1] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = -x * s + g * c;
                h = y * s;
                y = y * c;
                for (j = 0; j < n; j++) {
                    x = v[j][i - 1];
                    z = v[j][i];
                    v[j][i - 1] = x * c + z * s;
                    v[j][i] = -x * s + z * c
                }
                z = pythag(f, h);
                q[i - 1] = z;
                c = f / z;
                s = h / z;
                f = c * g + s * y;
                x = -s * g + c * y;
                for (j = 0; j < m; j++) {
                    y = u[j][i - 1];
                    z = u[j][i];
                    u[j][i - 1] = y * c + z * s;
                    u[j][i] = -y * s + z * c
                }
            }
            e[l] = 0;
            e[k] = f;
            q[k] = x
        }
    }
    for (i = 0; i < q.length; i++)
        if (q[i] < prec) q[i] = 0;
    for (i = 0; i < n; i++) {
        for (j = i - 1; j >= 0; j--) {
            if (q[j] < q[i]) {
                c = q[j];
                q[j] = q[i];
                q[i] = c;
                for (k = 0; k < u.length; k++) {
                    temp = u[k][i];
                    u[k][i] = u[k][j];
                    u[k][j] = temp
                }
                for (k = 0; k < v.length; k++) {
                    temp = v[k][i];
                    v[k][i] = v[k][j];
                    v[k][j] = temp
                }
                i = j
            }
        }
    }
    return {
        U: u,
        S: q,
        V: v
    }
};
'use strict';
angular.module('farmdoc').config(['$stateProvider', '$urlRouterProvider', function($stateProvider, $urlRouterProvider) {
    $urlRouterProvider.otherwise('/');
    $stateProvider.state('index', {
        url: '/',
        templateUrl: 'templates/premcalc.html'
    }).state('simulator', {
        url: '/simulator',
        templateUrl: 'templates/simulator.html'
    }).state('evaluator', {
        url: '/evaluator',
        templateUrl: 'templates/simulator.html'
    }).state('optionquery', {
        url: '/optionquery',
        templateUrl: 'templates/optionquery.html'
    }).state('arctool', {
        url: '/arctool',
        templateUrl: 'templates/arccalc.html'
    });
}]);
angular.module('farmdoc').controller('AlertsCtrl', ['$scope', AlertsCtrl]);

function AlertsCtrl($scope) {
    $scope.alerts = [];
    $scope.addAlert = function() {
        $scope.alerts.push({
            msg: 'Another alert!'
        });
    };
    $scope.closeAlert = function(index) {
        $scope.alerts.splice(index, 1);
    };
}
angular.module("farmdoc").controller('ARCCtrl', ['$scope', '$cookieStore', '$resource', '$location', '$anchorScroll', '$modal', 'config', 'nav', ARCCtrl]);

function ARCCtrl($scope, $cookieStore, $resource, $location, $anchorScroll, $modal, config, nav) {
    $scope.selection = {};
    nav.setPageName("ARC-CO Tool");
    initOptions();

    function initOptions() {
        $resource(config.serverRoot + '/data/init/simulator').get({}, function(response) {
            $scope.states = response.states;
            $scope.crops = response.crops;
            $scope.counties = response.ilcounty;
            $scope.selection.crop = "Corn";
            $scope.PLCyield = 0;
            $scope.selection.state = "Illinois";
            $scope.selection.county = "LaSalle";
            $scope.inputsChanged('Illinois', 'LaSalle', 'Corn', 'All');
        });
    }
    $scope.r2 = function(val) {
        return parseFloat(val).toFixed(2);
    }
    $scope.changedState = function() {
        $resource(config.serverRoot + '/data/countiesbystatename/:state').get({
            state: $scope.selection.state
        }, function(response) {
            $scope.counties = response.counties;
            $scope.selection.county = $scope.counties[0].ctyName;
        });
    };
    $scope.inputsChanged = function(state, county, crop, type) {
        state = $scope.selection.state;
        county = $scope.selection.county;
        crop = $scope.selection.crop;
        type = "All";
        $resource(config.serverRoot + "/arctool/lookup/" + state + "/" + county + "/" + crop + "/" + type).get({}, function(response) {
            $scope.results = response.results;
            $scope.yieldy2015 = $scope.results.arc_lookup.y2015;
            $scope.yieldy2016 = $scope.results.arc_lookup.y2016;
            $scope.yieldy2017 = $scope.results.arc_lookup.y2017;
            $scope.yieldy2018 = $scope.results.arc_lookup.y2018;
            $scope.pricey2015 = parseFloat(response.results.arc_price_year.y2015).toFixed(2);
            $scope.pricey2016 = parseFloat(response.results.arc_price_year.y2016).toFixed(2);
            $scope.pricey2017 = parseFloat(response.results.arc_price_year.y2017).toFixed(2);
            $scope.pricey2018 = parseFloat(response.results.arc_price_year.y2018).toFixed(2);
            $scope.PLCyield = Math.round(((parseInt($scope.results.arc_lookup.y2009) + parseInt($scope.results.arc_lookup.y2010) + parseInt($scope.results.arc_lookup.y2011) + parseInt($scope.results.arc_lookup.y2012)) / 4) * .85);
        });
    };
    $scope.calculateResults = function(number1, number2, number3, number4, number5) {
        var total = parseInt(number1) + parseInt(number2) + parseInt(number3) + parseInt(number4) + parseInt(number5);
        var min = Math.min(number1, number2, number3, number4, number5);
        var max = Math.max(number1, number2, number3, number4, number5);
        var average = (total - min - max) / 3;
        return Math.round(average);
    }
    $scope.calculateBenchPrice = function(arg1, arg2, arg3, arg4, arg5, arg6) {
        var price1 = parseFloat(arg6);
        var number1 = Math.max(parseFloat(arg1), price1);
        var number2 = Math.max(parseFloat(arg2), price1);
        var number3 = Math.max(parseFloat(arg3), price1);
        var number4 = Math.max(parseFloat(arg4), price1);
        var number5 = Math.max(parseFloat(arg5), price1);
        var total = number1 + number2 + number3 + number4 + number5;
        var min = Math.min(number1, number2, number3, number4, number5);
        var max = Math.max(number1, number2, number3, number4, number5);
        var average = (total - min - max) / 3;
        return average.toFixed(2);
    }
    $scope.calculateARCGuarantee = function(yield, price, display) {
        if (display == null) {
            display = 0;
        }
        var guarantee = .86 * yield* price;
        if (display == 1) {
            return Math.round(guarantee);
        } else {
            return guarantee.toFixed(2);
        }
    }
    $scope.calculateARCPrice = function(benchyield, benchprice, arcguarantee, countyyield, myaprice, price2, timesnumber) {
        if (timesnumber == null) {
            timesnumber = .85;
        }
        var maxnumbermatch = Math.max(parseFloat(myaprice), parseFloat(price2));
        var firstnumber = benchyield * benchprice * .1;
        var roundyieldresult = parseFloat(countyyield) * maxnumbermatch;
        var secondnumber = arcguarantee - roundyieldresult;
        var minnumbermatch = Math.min(firstnumber, secondnumber);
        var arcprice = Math.max(0, minnumbermatch) * timesnumber;
        return arcprice.toFixed(2);
    }
    $scope.calculatePLCPayment = function(plcyield, price1, MYAPrice, price2, timesnumber) {
        if (timesnumber == null) {
            timesnumber = .85;
        }
        var price2max = Math.max(MYAPrice, parseFloat(price2));
        var maxnumber1 = Math.max(0, parseFloat(price1) - price2max);
        var maxnumber2 = Math.max(0, maxnumber1.toFixed(4));
        var plcpayment = plcyield * timesnumber * maxnumber2;
        return Math.round(plcpayment);
    }
    $scope.calculateARCaverage = function(arc1, arc2, arc3, arc4, arc5) {
        var returnnumber = (parseFloat(arc1) + parseFloat(arc2) + parseFloat(arc3) + parseFloat(arc4) + parseFloat(arc5)) / 5;
        return returnnumber.toFixed(2);
    }
    $scope.calculatePLCaverage = function(arc1, arc2, arc3, arc4, arc5) {
        var returnnumber = (parseFloat(arc1) + parseFloat(arc2) + parseFloat(arc3) + parseFloat(arc4) + parseFloat(arc5)) / 5;
        return Math.round(returnnumber);
    }
    $scope.showast = function(a19) {
        if (a19 == 1) {
            return "*";
        } else {
            return "";
        }
    }
    $scope.calculate = function() {}

    function calculate() {}
}
angular.module('farmdoc').controller('MasterCtrl', ['$scope', '$cookieStore', 'nav', MasterCtrl]);

function MasterCtrl($scope, $cookieStore, nav) {
    var mobileView = 992;
    $scope.getWidth = function() {
        return window.innerWidth;
    };
    $scope.$watch($scope.getWidth, function(newValue, oldValue) {
        if (newValue >= mobileView) {
            if (angular.isDefined($cookieStore.get('toggle'))) {
                $scope.toggle = !$cookieStore.get('toggle') ? false : true;
            } else {
                $scope.toggle = true;
            }
        } else {
            $scope.toggle = false;
        }
    });
    $scope.toggleSidebar = function() {
        $scope.toggle = !$scope.toggle;
        $cookieStore.put('toggle', $scope.toggle);
    };
    window.onresize = function() {
        $scope.$apply();
    };
    $scope.getPageName = function() {
        return nav.getPageName();
    }
}
angular.module('farmdoc').controller('OptionQueryCtrl', ['$scope', '$cookieStore', '$resource', '$location', '$anchorScroll', '$modal', 'config', 'nav', OptionQueryCtrl]);

function OptionQueryCtrl($scope, $cookieStore, $resource, $location, $anchorScroll, $modal, config, nav) {
    nav.setPageName("Option Query");
    $scope.selection = {
        crop: 'corn',
        monthCode: 'Z',
        year: '17'
    };
    $scope.cropOptions = [{
        value: 'corn',
        displayName: 'Corn'
    }, {
        value: 'soybeans',
        displayName: 'Soybeans'
    }]
    $scope.inputs = {
        resultPercentiles: [5, 15, 25, 35, 45, 50, 55, 65, 75, 85, 95],
        priceOfInterest: 3.86
    };
    $scope.yearOptions = [{
        value: 17,
        displayName: 2017
    }, {
        value: 18,
        displayName: 2018
    }]
    $scope.monthCodes = {
        'corn': [{
            code: 'H',
            label: 'Mar'
        }, {
            code: 'K',
            label: 'May'
        }, {
            code: 'N',
            label: 'Jul'
        }, {
            code: 'U',
            label: 'Sep'
        }, {
            code: 'Z',
            label: 'Dec'
        }],
        'soybeans': [{
            code: 'F',
            label: 'Jan'
        }, {
            code: 'H',
            label: 'Mar'
        }, {
            code: 'K',
            label: 'May'
        }, {
            code: 'N',
            label: 'Jul'
        }, {
            code: 'Q',
            label: 'Aug'
        }, {
            code: 'U',
            label: 'Sep'
        }, {
            code: 'X',
            label: 'Nov'
        }]
    };
    $scope.newmonthCodes = function(crop, year) {
        var numbermonths = ["F", null, "H", null, "K", null, "N", "Q", "U", null, "X", "Z"];
        var dateObj = new Date();
        var month = dateObj.getUTCMonth() + 1;
        var day = dateObj.getUTCDate();
        var fullyear = dateObj.getUTCFullYear();
        var twodateyear = fullyear.toString().substring(2);
        if (year == twodateyear) {
            var mcodes = $scope.monthCodes[crop];
            var returncodes = [];
            for (i = month; i < numbermonths.length; i++) {
                var searchmonth = numbermonths[i];
                for (inc = 0; inc < mcodes.length; inc++) {
                    var item = mcodes[inc];
                    if (item.code == searchmonth) {
                        returncodes.push(item);
                    }
                }
            }
            return returncodes;
        } else {
            return $scope.monthCodes[crop];
        }
    }
    $scope.returnMonthLabel = function(cropCode) {
        if (cropCode == 'F') {
            return "Jan";
        }
        if (cropCode == 'H') {
            return "Mar";
        }
        if (cropCode == 'K') {
            return "May";
        }
        if (cropCode == 'N') {
            return "Jul";
        }
        if (cropCode == 'Q') {
            return "Aug";
        }
        if (cropCode == 'U') {
            return "Sep";
        }
        if (cropCode == 'X') {
            return "Nov";
        }
        if (cropCode == 'Z') {
            return "Dec";
        }
    }
    $scope.objectKeys = function(obj) {
        return Object.keys(obj);
    }
    $scope.status = {};

    function getPriceOfInterest(crop) {
        if (crop == 'soybeans') {
            return 8.85;
        }
        return 3.86;
    }
    $scope.date = new Date();
    $scope.changedCrop = function(crop, monthCode, year) {
        if (monthCode == null) {
            if (crop == 'corn') {
                $scope.selection.monthCode = 'Z';
                monthCode = 'Z';
            } else {
                $scope.selection.monthCode = 'X';
                monthCode = 'X';
            }
        }
        year = $scope.selection.year;
        if (year == null) {
            var dateObj = new Date();
            var fullyear = dateObj.getUTCFullYear();
            if (monthCode == 'Z' || monthCode == 'X') {
                fullyear++;
            }
            var twodateyear = fullyear.toString().substring(2);
            year = twodateyear;
        }
        $scope.status.error = null;
        $scope.status.loading = true;
        $scope.inputs.priceOfInterest = getPriceOfInterest(crop);
        $resource(config.serverRoot + "/optionquery/futuresoption/" + crop + "/" + monthCode + "/" + year).get({}, function(response) {
            $scope.results = response.results;
            var optimizationFunction = function(x) {
                return computeScenario(x[0], x[1], response.results.futuresData);
            };
            var solution = numeric.uncmin(optimizationFunction, [$scope.results.previousSolution.mu, $scope.results.previousSolution.sigma]).solution;
            $scope.inputs.priceOfInterest = $scope.results.price.toFixed(2);
            $scope.results.mu = solution[0];
            $scope.results.sigma = solution[1];
            var chartData = generateChartData($scope.results.sigma, $scope.results.mu);
            $scope.results.chartData = chartData;
            $scope.results.sigSquared = Math.pow($scope.results.sigma, 2);
            if ($scope.error) {
                $scope.status.error = $scope.error;
            } else {
                setTimeout(prepareCharts, 250);
            }
            $scope.regeneratePriceTableData();
            $resource(config.serverRoot + "/optionquery/solution").save({
                crop: crop,
                sigma: solution[0],
                mu: solution[1]
            });
            $scope.status.loading = false;
        });
    }
    $scope.changedCrop($scope.selection.crop, $scope.selection.monthCode, $scope.selection.year);
    $scope.regeneratePriceTableData = function() {
        var sigma = $scope.results.sigma;
        var mu = $scope.results.mu;
        var futuresPrice = $scope.results.price;
        var roundedToQuarter = (Math.round(futuresPrice * 4) / 4).toFixed(2);
        var current = roundedToQuarter - 1;
        $scope.results.priceTableData = [];
        for (var i = 0; i < 9; i++) {
            $scope.results.priceTableData.push({
                price: current,
                probability: getProbabilityForPrice(sigma, mu, current) * 100
            });
            current += .25;
        }
    }
    $scope.changedPriceofInterest = function() {
        $scope.regeneratePriceTableData();
        prepareCharts();
    }
    $scope.getProbabilityForPriceofInterest = function() {
        var sigma = $scope.results.sigma;
        var mu = $scope.results.mu;
        var poi = $scope.inputs.priceOfInterest;
        var probability = getProbabilityForPrice(sigma, mu, poi) * 100;
        return probability;
    }

    function getProbabilityForPrice(sigma, mu, price) {
        return normalDistribution((Math.log(price) - mu) / sigma);
    }

    function generateChartData(sigma, mu) {
        var ll = 0.0001;
        var ul = 0.999;
        var AB3 = Math.exp((normsInv(ll) * sigma + mu));
        var AB4 = Math.exp((normsInv(ul) * sigma + mu));
        var inc = (AB4 - AB3) / 599;
        var results = [];
        var current = AB3;
        for (i = 0; i < 600; i++) {
            var bigPyt = normalDistribution((Math.log(current) - mu) / sigma);
            var pi = Math.PI;
            var littlepyt = Math.exp((-0.5 * Math.pow(((Math.log(current) - mu) / sigma), 2))) / ((current * sigma * Math.pow(2 * pi, 0.5)));
            results.push({
                'price': current,
                'bigPyt': bigPyt,
                'littlepyt': littlepyt
            });
            current += inc;
        }
        return results;
    }
    $scope.computeBreakpoint = function(percentile) {
        percentile /= 100;
        return Math.exp((normsInv(percentile) * $scope.results.sigma + $scope.results.mu));
    };
    $scope.calcBigF = function(x) {
        return normalDistribution((Math.log(x) - $scope.results.mu) / $scope.results.sigma);
    };
    $scope.calcLittlef = function(x) {
        return Math.exp((-0.5 * Math.pow((Math.log(x) - $scope.results.mu) / $scope.results.sigma, 2))) / ((x * $scope.results.sigma * Math.pow(2 * Math.PI, 0.5)));
    };
    $scope.calcV = function(x) {
        return Math.exp(2 * $scope.results.mu + $scope.results.sigSquared) * (Math.exp($scope.results.sigSquared) - 1);
    };
    $scope.calcSig = function(x) {
        return Math.pow(Math.exp(2 * $scope.results.mu + $scope.results.sigSquared) * (Math.exp($scope.results.sigSquared) - 1), 0.5);
    };
    $scope.calcOptionBasedEx = function() {
        return Math.exp($scope.results.mu + ($scope.results.sigma * $scope.results.sigma) / 2);
    }
    $scope.calcIVRemaining = function(x) {
        return $scope.calcSig(x) / $scope.calcOptionBasedEx();
    };
    $scope.calcIVAnnualized = function(x) {
        return Math.pow(($scope.calcV(x) * 360 / $scope.results.dte), 0.5) / $scope.calcOptionBasedEx();
    }

    function prepareCharts() {
        prepareProbChart('Cumulative Probability of Prices at Expiration', 'Probability', 'bigPyt', 'chart1');
        prepareProbChart('Probability of Prices at Expiration', 'Relative Probability', 'littlepyt', 'chart2');
    }

    function prepareChartOptionValuesByStrike() {
        var rows = [];
        var dataTable = new google.visualization.DataTable();
        dataTable.addColumn('number', 'Strike Price');
        dataTable.addColumn('number', 'Call Premium');
        dataTable.addColumn('number', 'Put Premium');
        angular.forEach($scope.results.optionValuesByStrike, function(data) {
            var values = [parseFloat(data['strike']), parseFloat(data['put']), parseFloat(data['call'])];
            rows.push(values);
        });
        dataTable.addRows(rows);
        var moneyFormatter = new google.visualization.NumberFormat({
            fractionDigits: 2,
            prefix: '$'
        });
        moneyFormatter.format(dataTable, 0);
        var options = {
            title: 'Option Values by Strike',
            curveType: 'function',
            legend: {
                position: 'right',
                textStyle: {
                    fontSize: 14
                }
            },
            hAxis: {
                format: 'currency'
            },
            vAxis: {
                viewWindow: {
                    min: 0
                }
            },
            chartArea: {
                width: "50%",
                height: "80%"
            }
        };
        var chart = new google.visualization.LineChart(document.getElementById('chart3'));
        chart.draw(dataTable, options);
    }

    function prepareProbChart(title, yAxisLabel, dataColumn, chartDiv) {
        var rows = [];
        var poi = $scope.inputs.priceOfInterest;
        var dataTable = new google.visualization.DataTable();
        dataTable.addColumn('number', 'Price ($/bu)');
        dataTable.addColumn('number', yAxisLabel);
        dataTable.addColumn('number', 'POI line');
        angular.forEach($scope.results.chartData, function(data) {
            var poiLine = data[dataColumn];
            if (data['price'] > poi) {
                poiLine = 0;
            }
            var values = [parseFloat(data['price']), parseFloat(data[dataColumn]), parseFloat(poiLine)];
            rows.push(values);
        });
        dataTable.addRows(rows);
        var percentFormatter = new google.visualization.NumberFormat({
            fractionDigits: 2,
            suffix: '%'
        });
        var moneyFormatter = new google.visualization.NumberFormat({
            fractionDigits: 2,
            prefix: '$'
        });
        moneyFormatter.format(dataTable, 0);
        percentFormatter.format(dataTable, 1);
        var chartTitle = 'Probability';
        if (chartDiv == "chart2") {
            chartTitle = 'Relative Probability';
        }
        var options = {
            title: title,
            curveType: 'function',
            legend: {
                position: 'none',
                textStyle: {
                    fontSize: 14
                }
            },
            hAxis: {
                format: 'currency',
                slantedText: true,
                slantedTextAngle: 90,
                gridlines: {
                    count: 15
                }
            },
            vAxis: {
                title: chartTitle,
                viewWindow: {
                    min: 0
                },
                format: 'percent'
            },
            chartArea: {
                width: "80%",
                height: "68%"
            },
            colors: ['#C2D1F0', '#5b84d7']
        };
        var chart = new google.visualization.AreaChart(document.getElementById(chartDiv));
        chart.draw(dataTable, options);
    }

    function computeScenario(mu, sig, futuresData) {
        var cErrors = 0;
        var pErrors = 0;
        var below = 6;
        var above = 6;
        var step = 10;
        var nearest = Math.round((parseFloat(futuresData['price'])) / step) * step;
        var lower = nearest - (below * step);
        var upper = nearest + (above * step);
        var strike = lower;
        while (strike <= upper) {
            var rf = 0.03;
            var dte = futuresData['dte'];
            var key = "" + strike + "-0C";
            var callPrem = futuresData.data[key];
            var callDevSquared = 0;
            if (callPrem != null) {
                var btCalc = 1 / Math.exp(rf * dte / 360);
                var optionBasedEx = Math.exp(mu + (sig * sig / 2));
                var fk = normalDistribution((Math.log(strike / 100) - mu) / sig);
                var iexk = optionBasedEx * normalDistribution((Math.log(strike / 100) - mu - sig * sig) / sig);
                var callImpliedPrice = (optionBasedEx - iexk - strike / 100 * (1 - fk)) * btCalc;
                callDevSquared = Math.pow(callPrem / 100 - callImpliedPrice, 2);
            }
            var putDevSquared = 0;
            var putKey = "" + strike + "-0P";
            var putPrem = futuresData.data[putKey];
            if (putPrem != null) {
                var putPrice = ((strike / 100) * fk - iexk) * btCalc;
                putDevSquared = Math.pow(putPrem / 100 - putPrice, 2);
            }
            strike += step;
            cErrors += callDevSquared;
            pErrors += putDevSquared;
        }
        var totalErrors = cErrors + pErrors;
        return totalErrors;
    }

    function normalDistribution(z) {
        var c1 = 2.506628;
        var c2 = 0.3193815;
        var c3 = -0.3565638;
        var c4 = 1.7814779;
        var c5 = -1.821256;
        var c6 = 1.3302744;
        var w = 0;
        if (z > 0 || z == 0) {
            w = 1;
        } else {
            w = -1;
        }
        var y = 1 / (1 + 0.2316419 * w * z);
        return 0.5 + w * (0.5 - (Math.exp(-z * z / 2) / c1) * (y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * c6))))));
    }

    function normsInv(p) {
        if (typeof p == 'string' || p instanceof String) {
            p = parseFloat(p);
        }
        var a1 = -39.6968302866538;
        var a2 = 220.946098424521;
        var a3 = -275.928510446969;
        var a4 = 138.357751867269;
        var a5 = -30.6647980661472;
        var a6 = 2.50662827745924;
        var b1 = -54.4760987982241;
        var b2 = 161.585836858041;
        var b3 = -155.698979859887;
        var b4 = 66.8013118877197;
        var b5 = -13.2806815528857;
        var c1 = -.00778489400243029;
        var c2 = -0.322396458041136;
        var c3 = -2.40075827716184;
        var c4 = -2.54973253934373;
        var c5 = 4.37466414146497;
        var c6 = 2.93816398269878;
        var d1 = .00778469570904146;
        var d2 = 0.32246712907004;
        var d3 = 2.445134137143;
        var d4 = 3.75440866190742;
        var p_low = 0.02425;
        var p_high = 1 - p_low;
        var q = 0.0;
        var r = 0.0;
        var NormSInv = 0;
        if (p < 0 || p > 1) {
            return null;
        } else if (p < p_low) {
            q = Math.sqrt(-2 * Math.log(p));
            NormSInv = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
        } else if (p <= p_high) {
            q = p - 0.5;
            r = q * q;
            NormSInv = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
        } else {
            q = Math.sqrt(-2 * Math.log(1 - p));
            NormSInv = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
        }
        return NormSInv;
    }
}
angular.module('farmdoc').controller('PremCalcCtrl', ['$scope', '$cookieStore', '$resource', '$location', '$anchorScroll', '$modal', 'config', 'nav', PremCalcCtrl]);

function PremCalcCtrl($scope, $cookieStore, $resource, $location, $anchorScroll, $modal, config, nav) {
    $scope.levels = [50, 55, 60, 65, 70, 75, 80, 85];
    $scope.selection = {};
    $scope.counties = {};
    nav.setPageName("Premium Calculator");
    $resource(config.serverRoot + '/data/init/premcalc').get({}, function(response) {
        $scope.states = response.states;
        $scope.counties = response.ilcounty;
        $scope.crops = response.crops;
        $scope.selection.state = "17";
        $scope.selection.county = "1";
        $scope.selection.crop = "41";
        $scope.asof = response.asof;
        recalcParams();
    });
    $scope.changedState = function() {
        $scope.counties = {};
        $resource(config.serverRoot + '/data/counties/:state').get({
            state: $scope.selection.state
        }, function(response) {
            $scope.counties = response.counties;
            recalcParams();
        });
    };
    $scope.changedCounty = function() {
        recalcParams();
    };
    $scope.changedCrop = function() {
        recalcParams();
    };
    $scope.isCalculateReady = function() {
        return ($scope.counties.length > 0);
    };
    $scope.cropCode = function() {
        var sel = $scope.selection;
        if (!sel.crop) {
            return ''
        }
        if (!sel.county) {
            return '';
        }
        return (parseInt(sel.state) * 1000 + parseInt(sel.county)) * 100 + parseInt(sel.crop);
    };
    $scope.fullCode = function() {
        var code = $scope.cropCode();
        code = code * 1000 + parseInt($scope.params.type);
        code = code * 1000 + parseInt($scope.params.practice);
        return code;
    };
    $scope.countyFullCode = function() {
        var code = $scope.cropCode();
        code = code * 1000 + parseInt($scope.params.grpType);
        code = code * 1000 + parseInt($scope.params.grpPractice);
        return code;
    }
    $scope.showTaCalculator = function() {
        var modalInstance = $modal.open({
            templateUrl: 'templates/taWorksheet.html',
            controller: 'TaWorksheetCtrl',
            scope: $scope,
            size: 'lg'
        });
        modalInstance.result.then(function(result) {
            $scope.params.TAYield = result;
            $scope.params.useTaAdjustment = 1;
            $scope.calculateIndividual();
            $scope.calculateCounty();
        });
    };
    $scope.calculateAll = function() {
        $scope.calculateIndividual();
        $scope.calculateCounty();
    }
    $scope.calculateIndividual = function() {
        $scope.calculating = true;
        $scope.results = null;
        var calcParams = {
            code: $scope.fullCode(),
            comboProjPrice: $scope.params.comboProjPrice,
            acres: $scope.params.acres,
            comboVol: $scope.params.comboVol,
            hf: 0,
            preventedPlanting: $scope.params.preventedPlanting,
            riskMult: $scope.params.riskClass.riskMult,
            riskVal: $scope.params.riskClass.riskRate,
            TAYield: $scope.params.TAYield,
            useTaAdjustment: $scope.params.useTaAdjustment,
            aphYield: $scope.params.aphYield,
            rateYield: $scope.params.rateYield,
            riskClass: $scope.params.riskClass.riskType
        };
        $resource(config.serverRoot + '/compute/premiums').get(calcParams, function(response) {
            $scope.calculating = false;
            $scope.results = response.results;
            if ($scope.results.individual == null) {
                $scope.results.nodata = true;
            }
        });
    };
    $scope.calculateCounty = function() {
        $scope.calculatingCounty = true;
        $scope.resultsCounty = null;
        var calcParams = {
            code: $scope.countyFullCode(),
            vol: $scope.params.comboVol,
            bPrice: $scope.params.comboProjPrice
        };
        $resource(config.serverRoot + '/compute/premGrip').get(calcParams, function(response) {
            $scope.calculatingCounty = true;
            $scope.resultsCounty = response
            $scope.group = {
                percent: 95
            };
            if (response.grip.cov70 == 0 && response.gripHr.cov70 == 0 && response.grp.maxlib == 0) {
                $scope.resultsCounty.nodata = true;
            }
        });
    };

    function recalcParams() {
        var code = $scope.cropCode();
        if (code == '') {
            return;
        }
        $scope.params = {
            loading: true
        };
        $resource(config.serverRoot + '/compute/params/' + code + '/0/0').get({}, function(response) {
            $scope.params.loading = false;
            $scope.practices = response.practices;
            if ($scope.practices.length > 0) {
                $scope.params.practice = $scope.practices[0].practiceCode;
            }
            $scope.types = response.types;
            if ($scope.types.length > 0) {
                $scope.params.type = $scope.types[0].typeCode;
            }
            $scope.riskClasses = response.riskClasses;
            if ($scope.riskClasses.length > 0) {
                $scope.params.riskClass = $scope.riskClasses[0];
            }
            $scope.params.aphYield = response.aphYield;
            $scope.params.TAYield = response.TAYield;
            $scope.params.rateYield = response.rateYield;
            $scope.params.acres = response.acres;
            $scope.params.preventedPlanting = response.preventedPlanting;
            if (response.useTaAdjustment) {
                $scope.params.useTaAdjustment = 1;
            } else {
                $scope.params.useTaAdjustment = 0;
            }
            $scope.params.comboProjPrice = response.comboProjPrice;
            $scope.params.comboVol = response.comboVol;
            $scope.params.comboVolRevExc = response.comboVolRevExc;
            $scope.params.grpTypes = response.grpTypes;
            $scope.params.grpPractices = response.grpPractices;
            $scope.params.grpPrjPrice = response.grpPrjPrice;
            $scope.params.grpVol = response.grpVol;
            if ($scope.params.grpTypes.length > 0) {
                $scope.params.grpType = $scope.params.grpTypes[0].typeCode;
            }
            if ($scope.params.grpPractices.length > 0) {
                $scope.params.grpPractice = $scope.params.grpPractices[0].practiceCode;
            }
            $scope.calculateIndividual();
            $scope.calculateCounty();
        });
    }
}
angular.module('farmdoc').controller('SimulatorCtrl', ['$scope', '$cookieStore', '$resource', '$location', '$anchorScroll', '$modal', 'config', 'nav', SimulatorCtrl]);

function SimulatorCtrl($scope, $cookieStore, $resource, $location, $anchorScroll, $modal, config, nav) {
    $scope.selection = {};
    $scope.status = {
        loaded: false,
        calculating: false
    };
    nav.setPageName("Performance Evaluator");
    initOptions();

    function initOptions() {
        $resource(config.serverRoot + '/data/init/simulator').get({}, function(response) {
            $scope.states = response.states;
            $scope.crops = response.crops;
            $scope.counties = response.ilcounty;
            $scope.selection.crop = $scope.crops[0].cropCode;
            $scope.selection.unit = 'Basic';
            $scope.selection.acres = 320;
            $scope.selection.state = 17;
            $scope.selection.county = 1;
            $scope.asof = response.asof;
            calculate();
        });
    }
    $scope.r2 = function(val) {
        return parseFloat(val).toFixed(2);
    }
    $scope.changedState = function() {
        $scope.counties = {};
        $resource(config.serverRoot + '/data/counties/:state').get({
            state: $scope.selection.state
        }, function(response) {
            $scope.counties = response.counties;
        });
    };
    $scope.inputsChanged = function() {
        var grosstarget = "no";
        if ($scope.selection.state && $scope.selection.county) {
            calculate(grosstarget);
        }
    };
    $scope.grossTargetChanged = function() {
        var grosstarget = "yes";
        if ($scope.selection.state && $scope.selection.county) {
            calculate(grosstarget);
        }
    }
    $scope.cropCode = function() {
        var sel = $scope.selection;
        if (!sel.crop) {
            return ''
        }
        if (!sel.county) {
            return '';
        }
        return (parseInt(sel.state) * 1000 + parseInt(sel.county)) * 100 + parseInt(sel.crop);
    };
    $scope.scrollToGraph = function() {
        var chartTop = $('#curve_chart').offset().top;
        $(document.body).animate({
            scrollTop: chartTop
        });
        window.parent.parent.scrollTo(0, chartTop + 50);
    }

    function calculate(grosstarget) {
        $scope.status.loaded = false;
        $scope.status.calculating = true;
        if (grosstarget == "no") {
            $scope.selection.grossTarget = null;
        }
        var calcParams = {
            code: $scope.cropCode(),
            crop: $scope.selection.crop,
            enterpriseOrBasic: $scope.selection.unit,
            acres: $scope.selection.acres,
            grossTarget: $scope.selection.grossTarget
        };
        $resource(config.serverRoot + '/compute/simulator').get(calcParams, function(response) {
            $scope.status.calculating = false;
            $scope.status.loaded = true;
            $scope.selection.grossTarget = response.results.grossTarget;
            $scope.results = response.results;
            $scope.status.error = false;
            $scope.cfi = $scope.results.caseFarmInformation;
            if ($scope.results && $scope.results.estimatedPremiums && $scope.results.estimatedPremiums["50"].yp) {
                setTimeout(prepareChart, 250);
            } else {
                $scope.status.error = true;
            }
        });
    }

    function prepareChart() {
        var rows = [];
        var dataTable = new google.visualization.DataTable();
        dataTable.addColumn('number', 'Probability');
        dataTable.addColumn('number', 'No Ins');
        dataTable.addColumn('number', 'YP85');
        dataTable.addColumn('number', 'RP-HPE85');
        dataTable.addColumn('number', 'RP85');
        dataTable.addColumn('number', 'AYP90');
        dataTable.addColumn('number', 'ARPHPEP90');
        angular.forEach($scope.results.chartData, function(data, xAxis) {
            var values = [parseFloat(xAxis), parseFloat(data['No Ins']), parseFloat(data['YP85']), parseFloat(data['RP-HPE85']), parseFloat(data['RP85']), parseFloat(data['AYP90']), parseFloat(data['ARPHPE90'])];
            rows.push(values);
        });
        dataTable.addRows(rows);
        var percentFormatter = new google.visualization.NumberFormat({
            fractionDigits: 2,
            suffix: '%'
        });
        var moneyFormatter = new google.visualization.NumberFormat({
            fractionDigits: 2,
            prefix: '$'
        });
        moneyFormatter.format(dataTable, 0);
        percentFormatter.format(dataTable, 1);
        percentFormatter.format(dataTable, 2);
        percentFormatter.format(dataTable, 3);
        percentFormatter.format(dataTable, 4);
        percentFormatter.format(dataTable, 5);
        percentFormatter.format(dataTable, 6);
        var options = {
            title: 'Probabilities of Revenue with Insurance',
            curveType: 'function',
            legend: {
                position: 'right',
                textStyle: {
                    fontSize: 14
                }
            },
            hAxis: {
                format: 'currency'
            },
            vAxis: {
                viewWindow: {
                    min: 0
                },
                format: 'percent'
            },
            chartArea: {
                width: "50%",
                height: "80%"
            }
        };
        var chart = new google.visualization.LineChart(document.getElementById('curve_chart'));
        chart.draw(dataTable, options);
    }

    function drawChart() {
        var data = google.visualization.arrayToDataTable([
            ['Probability', 'No Ins.', 'YP85', 'RP-HPE85', 'RP85', 'AYP90', 'ARP90', 'ARPHE90'],
            [200, 0, 0, 0, 0, 0, 0, 0],
            [300, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            [400, 0.04, 0.04, 0.03, 0.015, 0.115, 0.115, 0.011],
            [500, 0.07, 0.06, 0.055, 0.20, 0.2, 0.133, 0.233],
            [600, 0.35, 0.35, 0.35, 0.25, 0.3, 0.233, 0.3]
        ]);
        var options = {
            title: 'Probabilities of Revenue with Insurance',
            curveType: 'function',
            legend: {
                position: 'bottom',
                textStyle: {
                    fontSize: 14
                }
            },
            hAxis: {
                format: 'currency',
            },
            vAxis: {
                viewWindow: {
                    min: 0
                }
            }
        };
        var chart = new google.visualization.LineChart(document.getElementById('curve_chart'));
        chart.draw(data, options);
    }
}
angular.module('farmdoc').controller('TaWorksheetCtrl', ['$scope', '$resource', '$modal', '$modalInstance', 'config', 'nav', TaWorksheetCtrl]);

function TaWorksheetCtrl($scope, $resource, $modal, $modalInstance, config, nav) {
    $scope.year = "2015";
    $scope.loading = {
        ta_init: true
    };
    $scope.yearData = [];
    $resource(config.serverRoot + '/data/yields/' + $scope.fullCode()).get({}, function(response) {
        var rawYieldData = response.yields;
        var currentYear = new Date().getFullYear();
        angular.forEach(rawYieldData, function(value) {
            if (value.year > currentYear - 12) {
                $scope.yearData.push({
                    year: value.year,
                    yield: value.yield,
                    yield_type: 'A'
                });
            }
        });
        $scope.loading = {};
        $scope.update();
    });
    $scope.update = function() {
        $scope.loading.ta_calc = true;
        $resource(config.serverRoot + "/compute/tayield/").get({
            code: $scope.fullCode(),
            yearData: angular.toJson($scope.yearData)
        }, function(response) {
            $scope.result = response;
            angular.forEach($scope.yearData, function(value) {
                var yearInfo = response.adjusted_data['year_' + value.year];
                value.yearly_adjustment = yearInfo.yearly_adjustment;
                value.ta_yield = yearInfo.ta_yield;
            });
            $scope.loading.ta_calc = false;
            if ($scope.result.ta_yield == 0 || $scope.result.ta_yield == null) {
                $scope.result.no_data = 1;
            }
        });
    }
    $scope.accept = function() {
        $scope.$close($scope.result.ta_yield);
    }
}
angular.module('farmdoc').directive('formOnChange', formOnChange);

function formOnChange($parse) {
    return {
        require: "form",
        link: function(scope, element, attrs) {
            var cb = $parse(attrs.formOnChange);
            element.on("change", function() {
                cb(scope);
            });
        }
    };
};
angular.module('farmdoc').directive('rdLoading', rdLoading);

function rdLoading() {
    var directive = {
        restrict: 'AE',
        template: '<div class="loading"><div class="double-bounce1"></div><div class="double-bounce2"></div></div>'
    };
    return directive;
};
angular.module('farmdoc').directive('simulatorBox', simulatorBox);

function simulatorBox() {
    var directive = {
        transclude: false,
        scope: {
            var: '=',
            prefix: '@',
            suffix: '@'
        },
        templateUrl: 'templates/simulatorBox.html',
        restrict: 'EA',
        link: function(scope, element, attrs) {
            scope.r2 = function(val) {
                if (!scope.prefix) {
                    scope.prefix = "";
                }
                if (!scope.suffix) {
                    scope.suffix = "";
                }
                if (parseFloat(val) == 0) {
                    return "NA";
                }
                if (scope.suffix == '%') {
                    return scope.prefix + "" + parseFloat(val * 100).toFixed(2) + scope.suffix;
                }
                return scope.prefix + "" + parseFloat(val).toFixed(2) + scope.suffix;
            }
        }
    };
    return directive;
}
angular.module('farmdoc').directive('simulatorBoxComplement', simulatorBoxComplement);

function simulatorBoxComplement() {
    var directive = {
        transclude: false,
        scope: {
            var: '=',
            prefix: '@',
            suffix: '@'
        },
        templateUrl: 'templates/simulatorBox.html',
        restrict: 'EA',
        link: function(scope, element, attrs) {
            scope.r2 = function(val) {
                if (!scope.prefix) {
                    scope.prefix = "";
                }
                if (!scope.suffix) {
                    scope.suffix = "";
                }
                if (parseFloat(val) == 0) {
                    return "NA";
                }
                if (scope.suffix == '%') {
                    return scope.prefix + "" + parseFloat((1 - val) * 100).toFixed(2) + scope.suffix;
                }
                return scope.prefix + "" + parseFloat(1 - val).toFixed(2) + scope.suffix;
            }
        }
    };
    return directive;
}
angular.module('farmdoc').directive('rdWidgetBody', rdWidgetBody);

function rdWidgetBody() {
    var directive = {
        requires: '^rdWidget',
        scope: {
            loading: '@?',
            classes: '@?'
        },
        transclude: true,
        template: '<div class="widget-body" ng-class="classes"><rd-loading ng-show="loading"></rd-loading><div ng-hide="loading" class="widget-content" ng-transclude></div></div>',
        restrict: 'E'
    };
    return directive;
};
angular.module('farmdoc').directive('rdWidgetFooter', rdWidgetFooter);

function rdWidgetFooter() {
    var directive = {
        requires: '^rdWidget',
        transclude: true,
        template: '<div class="widget-footer" ng-transclude></div>',
        restrict: 'E'
    };
    return directive;
};
angular.module('farmdoc').directive('rdWidgetHeader', rdWidgetTitle);

function rdWidgetTitle() {
    var directive = {
        requires: '^rdWidget',
        scope: {
            title: '@',
            icon: '@'
        },
        transclude: true,
        template: '<div class="widget-header"><div class="row"><div class="pull-left"><i class="fa" ng-class="icon"></i> {{title}} </div><div class="pull-right col-xs-6 col-sm-4" ng-transclude></div></div></div>',
        restrict: 'E'
    };
    return directive;
};
angular.module('farmdoc').directive('rdWidget', rdWidget);

function rdWidget() {
    var directive = {
        transclude: true,
        template: '<div class="widget" ng-transclude></div>',
        restrict: 'EA'
    };
    return directive;

    function link(scope, element, attrs) {}
};