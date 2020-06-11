// Extend the Array class
Array.prototype.max = function () {
    return Math.max.apply(null, this);
};
Array.prototype.min = function () {
    return Math.min.apply(null, this);
};
Array.prototype.mean = function () {
    var i, sum;
    for (i = 0, sum = 0; i < this.length; i++)
        sum += this[i];
    return sum / this.length;
};
Array.prototype.pip = function (x, y) {
    var i, j, c = false;
    for (i = 0, j = this.length - 1; i < this.length; j = i++) {
        if (((this[i][1] > y) != (this[j][1] > y)) &&
            (x < (this[j][0] - this[i][0]) * (y - this[i][1]) / (this[j][1] - this[i][1]) + this[i][0])) {
            c = !c;
        }
    }
    return c;
}

var kriging = function () {
    var kriging = {};

    //Matrix algebra
    var kriging_matrix_diag = function (c, n) {
        var Z = new Array(n), z, i, j;
        for (i = 0; i < n; i++) {
            z = new Array(n);
            for (j = 0; j < n; j++) {
                if (i == j) {
                    z[j] = c;
                } else {
                    z[j] = 0;
                }
            }
            Z[i] = z;
        }
        return Z;
    }

    var kriging_matrix_transpose = function (X, n, m) {
        var Z = new Array(m), z, i, j;
        for (i = 0; i < m; i++) {
            z = new Array(n);
            for (j = 0; j < n; j++) {
                z[j] = X[j][i];
            }
            Z[i] = z;
        }
        return Z;
    }

    var kriging_matrix_add = function (X, Y, n, m) {
        var Z = new Array(n), x, y, z, i, j;
        for (i = 0; i < n; i++) {
            x = X[i];
            y = Y[i];
            z = new Array(m);
            for (j = 0; j < m; j++) {
                z[j] = x[j] + y[j];
            }
            Z[i] = z;
        }
        return Z;
    }

    // Naive matrix multiplication
    var kriging_matrix_multiply = function (X, Y, n, m, p) {
        var Z = new Array(n), x, z, i, j, k;
        for (i = 0; i < n; i++) {
            x = X[i];
            z = new Array(p);
            for (j = 0; j < p; j++) {
                z[j] = 0;
                for (k = 0; k < m; k++) {
                    z[j] += x[k] * Y[k][j];
                }
            }
            Z[i] = z;
        }
        return Z;
    }

    // Cholesky decomposition
    var kriging_matrix_chol = function (X, n) {
        var p = new Array(n), x, i, j, k;
        for (i = 0; i < n; i++) {
            p[i] = X[i][i];
        }
        for (i = 0; i < n; i++) {
            x = X[i];
            for (j = 0; j < i; j++) {
                p[i] -= x[j] * x[j];
            }
            if (p[i] <= 0) {
                return false;
            }
            p[i] = Math.sqrt(p[i]);
            for (j = i + 1; j < n; j++) {
                for (k = 0; k < i; k++) {
                    X[j][i] -= X[j][k] * x[k];
                }
                X[j][i] /= p[i];
            }
        }
        for (i = 0; i < n; i++) {
            X[i][i] = p[i];
        }
        return true;
    }

    // Inversion of cholesky decomposition
    var kriging_matrix_chol2inv = function (X, n) {
        var sum, x, i, j, k;
        for (i = 0; i < n; i++) {
            X[i][i] = 1 / X[i][i];
            for (j = i + 1; j < n; j++) {
                sum = 0;
                for (k = i; k < j; k++) {
                    sum -= X[j][k] * X[k][i];
                }
                X[j][i] = sum / X[j][j];
            }
        }
        for (i = 0; i < n; i++) {
            x = X[i];
            for (j = i + 1; j < n; j++) {
                x[j] = 0;
            }
        }
        for (i = 0; i < n; i++) {
            x = X[i];
            x[i] *= x[i];
            for (j = i + 1; j < n; j++) {
                x[i] += X[j][i] * X[j][i];
            }
            for (j = i + 1; j < n; j++) {
                for (k = j; k < n; k++) {
                    x[j] += X[k][i] * X[k][j];
                }
            }
        }
        for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
                X[i][j] = X[j][i];
            }
        }
    }

    // Inversion via gauss-jordan elimination
    var kriging_matrix_solve = function (X, n) {
        var Y = kriging_matrix_diag(1, n);
        var indxc = new Array(n);
        var indxr = new Array(n);
        var ipiv = new Array(n);
        var i, j, k;
        for (i = 0; i < n; i++) {
            ipiv[i] = 0;
        }
        var irow, icol, temp, big, dum, pivinv;
        for (i = 0; i < n; i++) {
            big = 0;
            for (j = 0; j < n; j++) {
                if (ipiv[j] != 1) {
                    for (k = 0; k < n; k++) {
                        if (ipiv[k] == 0) {
                            if (Math.abs(X[j][k]) >= big) {
                                big = Math.abs(X[j][k]);
                                irow = j;
                                icol = k;
                            }
                        }
                    }
                }
            }
            ++(ipiv[icol]);
            if (irow != icol) {
                for (j = 0; j < n; j++) {
                    temp = X[irow][j];
                    X[irow][j] = X[icol][j];
                    X[icol][j] = temp;
                    temp = Y[irow][j];
                    Y[irow][j] = Y[icol][j];
                    Y[icol][j] = temp;
                }
            }
            indxr[i] = irow;
            indxc[i] = icol;
            if (X[icol][icol] == 0) {
                return false; // Singular
            }
            pivinv = 1 / X[icol][icol];
            X[icol][icol] = 1;
            for (j = 0; j < n; j++) {
                X[icol][j] *= pivinv;
                Y[icol][j] *= pivinv;
            }
            for (j = 0; j < n; j++) {
                if (j != icol) {
                    dum = X[j][icol];
                    X[j][icol] = 0;
                    for (k = 0; k < n; k++) {
                        X[j][k] -= X[icol][k] * dum;
                        Y[j][k] -= Y[icol][k] * dum;
                    }
                }
            }
        }
        for (i = n - 1; i >= 0; i--) {
            if (indxr[i] != indxc[i]) {
                for (j = 0; j < n; j++) {
                    temp = X[j][indxr[i]];
                    X[j][indxr[i]] = X[j][indxc[i]];
                    X[j][indxc[i]] = temp;
                }
            }
        }
        return true;
    }

    // Variogram models
    var kriging_variogram_gaussian = function (h, nugget, range, sill, A) {
        return nugget + ((sill - nugget) / range) * (1.0 - Math.exp(-(1.0 / A) * Math.pow(h / range, 2)));
    }

    var kriging_variogram_exponential = function (h, nugget, range, sill, A) {
        return nugget + ((sill - nugget) / range) * (1.0 - Math.exp(-(1.0 / A) * (h / range)));
    };

    var kriging_variogram_spherical = function (h, nugget, range, sill, A) {
        return h > range ? nugget + (sill - nugget) / range : nugget + ((sill - nugget) / range) * (1.5 * (h / range) - 0.5 * Math.pow(h / range, 3));
    };

    // Train using gaussian processes with bayesian priors
    kriging.train = function (t, x, y, model, sigma2, alpha) {
        var variogram = {
            t: t,
            x: x,
            y: y,
            nugget: 0.0,
            range: 0.0,
            sill: 0.0,
            A: 1 / 3,
            n: 0
        };
        switch (model) {
            case "gaussian": {
                variogram.model = kriging_variogram_gaussian;
                break;
            }
            case "exponential": {
                variogram.model = kriging_variogram_exponential;
                break;
            }
            case "spherical": {
                variogram.model = kriging_variogram_spherical;
                break;
            }
            default: {
                variogram.model = kriging_variogram_exponential;
                break;
            }
        }

        // Lag distance/semivariance
        var n = t.length;
        var distance = new Array((n * n - n) / 2);
        var i, j, k;
        for (i = 0, k = 0; i < n; i++) {
            for (j = i + 1; j < n; j++, k++) {
                distance[k] = new Array(2);
                distance[k][0] = Math.sqrt(Math.pow(x[i] - x[j], 2) + Math.pow(y[i] - y[j], 2));
                distance[k][1] = Math.abs(t[i] - t[j]);
            }
        }
        distance.sort(function (a, b) { return a[0] - b[0]; });
        variogram.range = distance[distance.length - 1][0];

        // Bin lag distance
        var lags = distance.length > 30 ? 30 : distance.length;
        var tolerance = variogram.range / lags;
        var lag = new Array(lags);
        var semi = new Array(lags);
        var m = 0, count;
        if (lags < 30) {
            for (; m < lags; m++) {
                lag[m] = distance[m][0];
                semi[m] = distance[m][1];
            }
        } else {
            for (i = 0, j = 0; i < lags && j < distance.length; i++) {
                lag[i] = 0;
                semi[i] = 0;
                count = 0;
                while (distance[j][0] <= (i + 1) * tolerance) {
                    lag[m] += distance[j][0];
                    semi[m] += distance[j][1];
                    j++;
                    count++;
                    if (j == distance.length) {
                        break;
                    }
                }
                if (count > 0) {
                    lag[m] /= count;
                    semi[m] /= count;
                    m++;
                }
            }
        }
        if (m < 2) {
            return variogram; // Error: Not enough points
        }

        // Feature transformation
        variogram.range = lag[m - 1] - lag[0];
        var X = new Array(m);
        var Y = new Array(m);
        var A = variogram.A;
        for (i = 0; i < m; i++) {
            X[i] = new Array(2);
            X[i][0] = 1;
            switch (model) {
                case "gaussian": {
                    X[i][1] = 1.0 - Math.exp(-(1.0 / A) * math.pow(lag[i] / variogram.range, 2));
                    break;
                }
                case "exponential": {
                    X[i][1] = 1.0 - Math.exp(-(1.0 / A) * lag[i] / variogram.range);
                    break;
                }
                case "spherical": {
                    X[i][1] = 1.5 * (lag[i] / variogram.range) - 0.5 * Math.pow(lag[i] / variogram.range, 3);
                    break;
                }
                default: {
                    X[i][1] = 1.0 - Math.exp(-(1.0 / A) * lag[i] / variogram.range);
                    break;
                }
            }
            Y[i] = [semi[i]];
        }

        // Least squares
        var Xt = kriging_matrix_transpose(X, m, 2);
        var Z = kriging_matrix_multiply(Xt, X, 2, m, 2);
        Z = kriging_matrix_add(Z, kriging_matrix_diag(1 / alpha, 2), 2, 2);
        var cloneZ = new Array(Z.length);
        for (i = 0; i < Z.length; i++) {
            cloneZ[i] = Z[i].slice(0);
        }
        if (kriging_matrix_chol(Z, 2)) {
            kriging_matrix_chol2inv(Z, 2);
        } else {
            kriging_matrix_solve(cloneZ, 2);
            Z = cloneZ;
        }
        var W = kriging_matrix_multiply(kriging_matrix_multiply(Z, Xt, 2, 2, m), Y, 2, m, 1);

        // Variogram parameters
        variogram.nugget = W[0][0];
        variogram.sill = W[1][0] * variogram.range + variogram.nugget;
        variogram.n = n;

        // Gram matrix with prior
        var K = new Array(n);
        for (i = 0; i < n; i++) {
            K[i] = new Array(n);
            for (j = 0; j < i; j++) {
                K[i][j] = variogram.model(Math.sqrt(Math.pow(x[i] - x[j], 2) + Math.pow(y[i] - y[j], 2)), variogram.nugget, variogram.range, variogram.sill, variogram.A);
                K[j][i] = K[i][j];
            }
            K[i][i] = variogram.model(0, variogram.nugget, variogram.range, variogram.sill, variogram.A);
        }

        // Inverse penalized Gram matrix projected to target vector
        var C = kriging_matrix_add(K, kriging_matrix_diag(sigma2, n), n, n);
        var cloneC = new Array(n);
        for (i = 0; i < n; i++) {
            cloneC[i] = C[i].slice(0);
        }
        if (kriging_matrix_chol(C, n)) {
            kriging_matrix_chol2inv(C, n);
        } else {
            kriging_matrix_solve(cloneC, n);
            C = cloneC;
        }

        // Copy unprojected inverted matrix as K
        var T = new Array(n);
        for (i = 0; i < n; i++) {
            T[i] = [t[i]];
        }
        var M = kriging_matrix_multiply(C, T, n, n, 1);
        variogram.K = C;
        variogram.M = M;

        return variogram;
    }

    kriging.predict = function (x, y, variogram) {
        var n = variogram.n;
        var k = new Array(n);
        for (var i = 0; i < n; i++) {
            k[i] = variogram.model(Math.sqrt(Math.pow(x - variogram.x[i], 2) + Math.pow(y - variogram.y[i], 2)), variogram.nugget, variogram.range, variogram.sill, variogram.A);
        }
        var K = [k];
        return kriging_matrix_multiply(K, variogram.M, 1, n, 1)[0][0];
    }

    // Gridded matrices
    kriging.grid = function (polygons, variogram, width) {
        var i, j, k, n = polygons.length;
        if (n == 0) return;

        // Boundaries of polygons space
        var xlim = [polygons[0][0][0], polygons[0][0][0]];
        var ylim = [polygons[0][0][1], polygons[0][0][1]];
        for (i = 0; i < n; i++) // Polygons
            for (j = 0; j < polygons[i].length; j++) { // Vertices
                if (polygons[i][j][0] < xlim[0])
                    xlim[0] = polygons[i][j][0];
                if (polygons[i][j][0] > xlim[1])
                    xlim[1] = polygons[i][j][0];
                if (polygons[i][j][1] < ylim[0])
                    ylim[0] = polygons[i][j][1];
                if (polygons[i][j][1] > ylim[1])
                    ylim[1] = polygons[i][j][1];
            }

        // Alloc for O(n^2) space
        var xtarget, ytarget;
        var a = Array(2), b = Array(2);
        var lxlim = Array(2); // Local dimensions
        var lylim = Array(2); // Local dimensions
        var x = Math.ceil((xlim[1] - xlim[0]) / width);
        var y = Math.ceil((ylim[1] - ylim[0]) / width);

        var A = Array(x + 1);
        for (i = 0; i <= x; i++) A[i] = Array(y + 1);
        for (i = 0; i < n; i++) {
            // Range for polygons[i]
            lxlim[0] = polygons[i][0][0];
            lxlim[1] = lxlim[0];
            lylim[0] = polygons[i][0][1];
            lylim[1] = lylim[0];
            for (j = 1; j < polygons[i].length; j++) { // Vertices
                if (polygons[i][j][0] < lxlim[0])
                    lxlim[0] = polygons[i][j][0];
                if (polygons[i][j][0] > lxlim[1])
                    lxlim[1] = polygons[i][j][0];
                if (polygons[i][j][1] < lylim[0])
                    lylim[0] = polygons[i][j][1];
                if (polygons[i][j][1] > lylim[1])
                    lylim[1] = polygons[i][j][1];
            }

            // Loop through polygon subspace
            a[0] = Math.floor(((lxlim[0] - ((lxlim[0] - xlim[0]) % width)) - xlim[0]) / width);
            a[1] = Math.ceil(((lxlim[1] - ((lxlim[1] - xlim[1]) % width)) - xlim[0]) / width);
            b[0] = Math.floor(((lylim[0] - ((lylim[0] - ylim[0]) % width)) - ylim[0]) / width);
            b[1] = Math.ceil(((lylim[1] - ((lylim[1] - ylim[1]) % width)) - ylim[0]) / width);
            for (j = a[0]; j <= a[1]; j++)
                for (k = b[0]; k <= b[1]; k++) {
                    xtarget = xlim[0] + j * width;
                    ytarget = ylim[0] + k * width;
                    if (polygons[i].pip(xtarget, ytarget))
                        A[j][k] = kriging.predict(xtarget,
                            ytarget,
                            variogram);
                }
        }
        A.xlim = xlim;
        A.ylim = ylim;
        A.zlim = [variogram.t.min(), variogram.t.max()];
        A.width = width;
        return A;
    };

    // Plotting on the DOM
    kriging.plot = function (canvas, grid, xlim, ylim, colors) {
        // Clear screen 
        var ctx = canvas.getContext("2d");
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // Starting boundaries
        var range = [xlim[1] - xlim[0], ylim[1] - ylim[0], grid.zlim[1] - grid.zlim[0]];
        var i, j, x, y, z;
        var n = grid.length;
        var m = grid[0].length;
        var wx = Math.ceil(grid.width * canvas.width / (xlim[1] - xlim[0]));
        var wy = Math.ceil(grid.width * canvas.height / (ylim[1] - ylim[0]));
        for (i = 0; i < n; i++)
            for (j = 0; j < m; j++) {
                if (grid[i][j] == undefined) continue;
                x = canvas.width * (i * grid.width + grid.xlim[0] - xlim[0]) / range[0];
                y = canvas.height * (1 - (j * grid.width + grid.ylim[0] - ylim[0]) / range[1]);
                //z = (grid[i][j] - grid.zlim[0]) / range[2];
                //if (z < 0.0) z = 0.0;
                //if (z > 1.0) z = 1.0;

                //ctx.fillStyle = colors[Math.floor((colors.length - 1) * z)];
                ctx.fillStyle = colors[Math.floor(grid[i][j] / 50)];
                ctx.fillRect(Math.round(x - wx / 2), Math.round(y - wy / 2), wx, wy);
            }

    };

    return kriging;
}();