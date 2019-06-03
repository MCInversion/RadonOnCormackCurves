#include "pch.h"
#include <iostream>
#include <Magick++.h>
#include <cmath>
#include <complex>
#include <cstdlib>

double polarToX(double r, double theta) {
	return r * cos(theta);
}

double polarToY(double r, double theta) {
	return r * sin(theta);
}

double XYtoR(double x, double y) {
	return sqrt(x * x + y * y);
}

double XYtoTheta(double x, double y) {
	return atan2(y, x);
}

int XtoPixelX(int w, double x) {
	return (int)round(x * w) + w / 2;
}

int YtoPixelY(int h, double y) {
	return (int)round(y * h) - y / 2;
}

bool isInsideImage(int x, int y, int h, int w) {
	return (x > 0 && x < w) && (y > 0 && y < h);
}

double f_gaussians(double x, double y) {
	return 2 * exp(-(x * x + y * y)) + exp(-((x - 3) * (x - 3) + (y - 3) * (y - 3)));
}

double f_gaussiand_polar(double r, double theta) {
	return f_gaussians(polarToX(r, theta), polarToY(r, theta));
}

void fft(double *x_in, std::complex<double> *x_out,	int N) {
	// Make copy of array and apply window
	for (int i = 0; i < N; i++) {
		x_out[i] = std::complex<double>(x_in[i], 0);
		x_out[i] *= 1; // Window
	}

	// Start recursion
	fft_rec(x_out, N);
}

void fft_rec(std::complex<double> *x, int N) {
	// Check if it is splitted enough
	if (N <= 1) {
		return;
	}

	// Split even and odd
	std::complex<double> *odd = new std::complex<double>[N / 2];
	std::complex<double> *even = new std::complex<double>[N / 2];
	for (int i = 0; i < N / 2; i++) {
		even[i] = x[i * 2];
		odd[i] = x[i * 2 + 1];
	}

	// Split on tasks
	fft_rec(even, N / 2);
	fft_rec(odd, N / 2);


	// Calculate DFT
	for (int k = 0; k < N / 2; k++) {
		std::complex<double> t = exp(std::complex<double>(0, -2 * M_PI * k / N)) * odd[k];
		x[k] = even[k] + t;
		x[N / 2 + k] = even[k] - t;
	}

	delete[] odd;
	delete[] even;
}

Magick::Image CormackTransform(Magick::Image *image, double sigma, bool alpha_transform = true) {
	const double phiMin = 0.; const double phiMax = 2 * M_PI;
	const double pMin = 0.; const double pMax = 1.;
	const double psiMax = (alpha_transform ? 0.99 : 1.) * M_PI / (2 * sigma);
	const double psiMin = -psiMax;
	double phi, p, dArcLength;
	double fvalR, fvalG, fvalB, fvalR_next, fvalG_next, fvalB_next;
	double psi, r, r_next, theta, theta_next, x, y, x_next, y_next, psiIntegSumR, psiIntegSumG, psiIntegSumB;
	int imgX, imgY, imgX_next, imgY_next;
	int imgP, imgPhi;
	Magick::ColorRGB px, px_next;

	const int w = image->columns();
	const int h = image->rows();

	// value arrays
	double* resValuesR = new double[w * h];
	double* resValuesG = new double[w * h];
	double* resValuesB = new double[w * h];

	double minR = DBL_MAX;
	double maxR = -DBL_MAX;
	double minG = minR;
	double maxG = maxR;
	double minB = minR;
	double maxB = maxR;

	const int integralSize = 400;
	double dPsi = (psiMax - psiMin) / integralSize;
	int percent = (h > 100 ? (int)round((double)h / 100.) : 1);

	std::string paramString = alpha_transform ? "alpha" : "beta";
	std::cout << "initializing Cormack " << paramString << " transform \\w " << paramString << " = " << sigma << std::endl;
	std::cout << "image : " << image->fileName() << ", " << h << "x" << w << std::endl;

	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			phi = (phiMax - phiMin) * i / h + phiMin;
			p = (pMax - pMin) * j / w + pMin;

			imgP = XtoPixelX((int)round((double)w / (pMax - pMin)), p);
			imgPhi = YtoPixelY((int)round((double)h / (phiMax - phiMin)), phi);

			psiIntegSumR = 0.; psiIntegSumG = 0.; psiIntegSumB = 0.;
			psi = psiMin;

			for (int k = 0; k < integralSize; k++) {
				psi += dPsi;

				r = (alpha_transform ? p / pow(cos(sigma * psi), 1. / sigma) : p * pow(cos(sigma * psi), 1. / sigma));
				r_next = (alpha_transform ? p / pow(cos(sigma * (psi + dPsi)), 1. / sigma) : p * pow(cos(sigma * (psi + dPsi)), 1. / sigma));
				theta = phi + psi;
				theta_next = theta + dPsi;

				x = polarToX(r, theta);
				y = polarToY(r, theta);

				x_next = polarToX(r_next, theta_next);
				y_next = polarToY(r_next, theta_next);

				imgX = XtoPixelX(w, x);
				imgY = YtoPixelY(h, y);

				imgX_next = XtoPixelX(w, x_next);
				imgY_next = YtoPixelY(h, y_next);

				if (!isInsideImage(imgX, imgY, h, w) || !isInsideImage(imgX_next, imgY_next, h, w)) {
					continue;
				}

				px = image->pixelColor(imgX, imgY);
				px_next = image->pixelColor(imgX_next, imgY_next);

				fvalR = px.red(); fvalG = px.green(); fvalB = px.blue();
				fvalR_next = px_next.red(); fvalG_next = px_next.green(); fvalB_next = px_next.blue();

				// fvalR = fvalG = fvalB = f_gaussians(x, y);

				dArcLength = (alpha_transform ? dPsi / pow(cos(sigma * psi), (1. + 1. / sigma)) : dPsi * pow(cos(sigma * psi), (1. / sigma - 1)));

				psiIntegSumR += p * 0.5 * (fvalR + fvalR_next) * dArcLength;
				psiIntegSumG += p * 0.5 * (fvalG + fvalG_next) * dArcLength;
				psiIntegSumB += p * 0.5 * (fvalB + fvalB_next) * dArcLength;
			}

			if (psiIntegSumR < minR) {
				minR = psiIntegSumR;
			}
			if (psiIntegSumR > maxR) {
				maxR = psiIntegSumR;
			}

			if (psiIntegSumG < minG) {
				minG = psiIntegSumG;
			}
			if (psiIntegSumG > maxG) {
				maxG = psiIntegSumG;
			}

			if (psiIntegSumB < minB) {
				minB = psiIntegSumB;
			}
			if (psiIntegSumB > maxB) {
				maxB = psiIntegSumB;
			}

			resValuesR[i * h + j] = psiIntegSumR;
			resValuesG[i * h + j] = psiIntegSumG;
			resValuesB[i * h + j] = psiIntegSumB;
		}


		if (i % percent == 0) {
			std::cout << "\rCormack Transform on " << image->fileName() << ": "
				<< (int)(((double)i / (double)h * 100) + 0.5) << " % complete ";
		}
	}

	std::cout << std::endl;
	std::cout << "RGB limits:" << std::endl;
	std::cout << "R: [" << minR << ", " << maxR << "]" << std::endl;
	std::cout << "G: [" << minG << ", " << maxG << "]" << std::endl;
	std::cout << "B: [" << minB << ", " << maxB << "]" << std::endl;

	// writing values into the result image
	Magick::Image result(Magick::Geometry(w, h), "black");
	result.type(Magick::TrueColorType);
	result.modifyImage();

	Magick::Pixels view(result);
	Magick::Quantum *pixels = view.get(0, 0, w, h);
	Magick::ColorRGB col;

	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			col = Magick::ColorRGB(
				resValuesR[i * h + j] / (maxR - minR) - minR,
				resValuesG[i * h + j] / (maxG - minG) - minG,
				resValuesB[i * h + j] / (maxB - minB) - minB
			);

			*pixels++ = col.quantumRed();
			*pixels++ = col.quantumGreen();
			*pixels++ = col.quantumBlue();
		}

		if (i % percent == 0) {
			std::cout << "\rWriting results : " << (int)(((double)i / (double)h * 100) + 0.5) << " % complete ";
		}
	}
	std::cout << std::endl;
	std::cout << "Cormack " << paramString << "-transform complete." << std::endl;


	view.sync();

	delete[] resValuesR;
	delete[] resValuesG;
	delete[] resValuesB;

	return result;
}

Magick::Image imageFFT(Magick::Image *image) {
	Magick::ColorRGB px;

	const int w = image->columns();
	const int h = image->rows();
	Magick::Image result(Magick::Geometry(w, h), "black");
	int percent = (h > 100 ? (int)round((double)h / 100.) : 1);

	result.type(Magick::TrueColorType);
	result.modifyImage();

	Magick::Pixels view(result);
	Magick::Quantum *pixels = view.get(0, 0, w, h);
	Magick::ColorRGB col;

	double *rowR = new double[w];
	double *rowG = new double[w];
	double *rowB = new double[w];

	std::complex<double> **outR = new std::complex<double>*[h];
	std::complex<double> **outG = new std::complex<double>*[h];
	std::complex<double> **outB = new std::complex<double>*[h];

	double minR = DBL_MAX;
	double maxR = -DBL_MAX;
	double minG = minR;
	double maxG = maxR;
	double minB = minR;
	double maxB = maxR;

	for (int i = 0; i < h; i++) {
		// fill 3 1D vectors
		outR[i] = new std::complex<double>[w];
		outG[i] = new std::complex<double>[w];
		outB[i] = new std::complex<double>[w];

		for (int j = 0; j < w; j++) {
			px = image->pixelColor(i, j);

			rowR[j] = px.red();
			rowG[j] = px.green();
			rowB[j] = px.blue();
		}
		fft(rowR, outR[i], w);
		fft(rowG, outG[i], w);
		fft(rowB, outB[i], w);

		for (int j = 0; j < w; j++) {
			if (outR[i][j].real() < minR) {
				minR = outR[i][j].real();
			}
			if (outR[i][j].real() > maxR) {
				maxR = outR[i][j].real();
			}

			if (outG[i][j].real() < minG) {
				minG = outG[i][j].real();
			}
			if (outG[i][j].real() > maxG) {
				maxG = outG[i][j].real();
			}

			if (outB[i][j].real() < minB) {
				minB = outB[i][j].real();
			}
			if (outB[i][j].real() > maxB) {
				maxB = outB[i][j].real();
			}
		}

		if (i % percent == 0) {
			std::cout << "\rFFT filling complex arrays : " << (int)(((double)i / (double)h * 100) + 0.5) << " % complete ";
		}
	}
	std::cout << std::endl;
	std::cout << "RGB limits:" << std::endl;
	std::cout << "R: [" << minR << ", " << maxR << "]" << std::endl;
	std::cout << "G: [" << minG << ", " << maxG << "]" << std::endl;
	std::cout << "B: [" << minB << ", " << maxB << "]" << std::endl;

	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			col = Magick::ColorRGB(
				outR[i][j].real() / (maxR - minR) - minR,
				outG[i][j].real() / (maxG - minG) - minG,
				outB[i][j].real() / (maxB - minB) - minB
			);

			*pixels++ = col.quantumRed();
			*pixels++ = col.quantumGreen();
			*pixels++ = col.quantumBlue();
		}

		if (i % percent == 0) {
			std::cout << "\rFFT translating to RGB : " << (int)(((double)i / (double)h * 100) + 0.5) << " % complete ";
		}

		delete[] outR[i]; delete[] outG[i]; delete[] outB[i];
	}
	view.sync();

	delete[] outR; delete[] outG; delete[] outB;
	delete[] rowR; delete[] rowG; delete[] rowB;

	std::cout << "image FFT complete." << std::endl;

	return result;
}

Magick::Image InverseCormackTransform(Magick::Image *image, double sigma, bool alpha_transform = true) {
	Magick::Image fourier_img = imageFFT(image);

	// fourier_img.inverseFourierTransform(*image, true);
	Magick::ColorRGB px, px_next;

	const int w = image->columns();
	const int h = image->rows();

	// value arrays
	double* resValuesR = new double[w * h];
	double* resValuesG = new double[w * h];
	double* resValuesB = new double[w * h];

	double minR = DBL_MAX;
	double maxR = -DBL_MAX;
	double minG = minR;
	double maxG = maxR;
	double minB = minR;
	double maxB = maxR;

	const int integralSize = 400;
	int percent = (w > 100 ? (int)round((double)w / 100.) : 1);
	int N = 5;

	double qIntegralSumR, qIntegralSumG, qIntegralSumB;
	double qIntegralSumR_next, qIntegralSumG_next, qIntegralSumB_next;
	double qIntegralR_dr, qIntegralG_dr, qIntegralB_dr;
	double fhatR, fhatG, fhatB, fhatR_next, fhatG_next, fhatB_next;
	double integrand, integrand_next;
	double cHarmonicSeriesSumR, cHarmonicSeriesSumG, cHarmonicSeriesSumB;
	double seriesRconst;
	double q, dq, theta;
	int imgQ, imgTheta, imgQ_next, imgTheta_next;
	double pMin = 0.01; double pMax = 1.;
	double thetaMin = 0.; double thetaMax = 2 * M_PI;
	double rMin = 0.01; double rMax = 1.;

	double p = pMin, dp = (pMax - pMin) / integralSize;
	double r = rMin;
	double dr = (rMax - rMin) / integralSize;
	double dTheta = (thetaMax - thetaMin) / integralSize;
	double qMin, qMax = 2.;

	std::string paramString = alpha_transform ? "alpha" : "beta";
	std::cout << "initializing Inverse Cormack " << paramString << "-transform \\w " << paramString << " = " << sigma << ", Fourier order: " << N << std::endl;
	std::cout << "image : " << h << "x" << w << std::endl;

	// for each (r, theta) we compute a circular harmonic sum of the q-integrals of Fourier image terms \hat{F}_l (q)
	for (int j = 0; j < w; j++) {
		qMin = pow(r, sigma);
		dq = (qMax - qMin) / integralSize;
		seriesRconst = -sigma * pow(r, sigma - 1) / M_PI;

		// angle cycle [0, 2 * Pi]
		theta = thetaMin;
		for (int i = 0; i < h; i++) {
			imgTheta = YtoPixelY(h, theta);

			// sum of N terms of circular harmonic decomposition of the Radon image
			cHarmonicSeriesSumR = 0.; cHarmonicSeriesSumG = 0.; cHarmonicSeriesSumB = 0.;
			for (int l = 0; l < N; l++) {
				// first integral for r
				qIntegralSumR = 0.; qIntegralSumG = 0.; qIntegralSumB = 0.;
				qIntegralSumR_next = 0.; qIntegralSumG_next = 0.; qIntegralSumB_next = 0.;

				q = qMin;
				for (int k = 0; k < integralSize; k++) {
					// \hat{F}_l (q) = \hat{f}_l(q^(1 / alpha))
					imgQ = XtoPixelX(w, pow(q, 1 / sigma));
					imgQ_next = XtoPixelX(w, pow(q + dq, 1 / sigma));

					if (!isInsideImage(imgQ, imgTheta, h, w) || !isInsideImage(imgQ_next, imgTheta, h, w)) {
						// std::cout << "\r (" << imgQ << ", " << imgTheta << ") or (" << imgQ_next << ", " << imgTheta << ") is outside image!";
						continue;
					} else {
						// std::cout << "\r (" << imgQ << ", " << imgTheta << ") or (" << imgQ_next << ", " << imgTheta << ") are inside";
					}

					px = fourier_img.pixelColor(imgQ, imgTheta);
					px_next = fourier_img.pixelColor(imgQ_next, imgTheta);

					fhatR = px.red();
					fhatG = px.green();
					fhatB = px.blue();

					fhatR_next = px_next.red();
					fhatG_next = px_next.green();
					fhatB_next = px_next.blue();

					integrand = cosh(l / sigma * acosh(q / pow(r, sigma)));
					integrand /= q * sqrt((q / pow(r, sigma)) * (q / pow(r, sigma)) - 1);

					integrand_next = cosh(l / sigma * acosh((q + dq) / pow(r, sigma)));
					integrand_next /= (q + dq)  * sqrt(((q + dq) / pow(r, sigma)) * ((q + dq) / pow(r, sigma)) - 1);

					qIntegralSumR += 0.5 * (integrand * fhatR + integrand_next * fhatR_next) * dq;
					qIntegralSumG += 0.5 * (integrand * fhatG + integrand_next * fhatG_next) * dq;
					qIntegralSumB += 0.5 * (integrand * fhatB + integrand_next * fhatB_next) * dq;

					q += dq;
				}

				qIntegralR_dr = qIntegralSumR;
				qIntegralG_dr = qIntegralSumG;
				qIntegralB_dr = qIntegralSumB;

				// second integral for r + dr
				qIntegralSumR = 0.; qIntegralSumG = 0.; qIntegralSumB = 0.;
				qIntegralSumR_next = 0.; qIntegralSumG_next = 0.; qIntegralSumB_next = 0.;

				q = qMin;
				for (int k = 0; k < integralSize; k++) {
					// \hat{F}_l (q) = \hat{f}_l(q^(1 / alpha))
					imgQ = XtoPixelX(w, pow(q, 1 / sigma));
					imgQ_next = XtoPixelX(w, pow(q + dq, 1 / sigma));

					if (!isInsideImage(imgQ, imgTheta, h, w) || !isInsideImage(imgQ_next, imgTheta, h, w)) {
						// std::cout << "\r (" << imgQ << ", " << imgTheta << ") or (" << imgQ_next << ", " << imgTheta << ") is outside image!";
						continue;
					}
					else {
						// std::cout << "\r (" << imgQ << ", " << imgTheta << ") or (" << imgQ_next << ", " << imgTheta << ") are inside";
					}

					// std::cout << "integral is being evaluated" << std::endl;

					px = fourier_img.pixelColor(imgQ, imgTheta);
					px_next = fourier_img.pixelColor(imgQ_next, imgTheta);

					fhatR = px.red();
					fhatG = px.green();
					fhatB = px.blue();

					fhatR_next = px_next.red();
					fhatG_next = px_next.green();
					fhatB_next = px_next.blue();

					integrand = cosh(l / sigma * acosh(q / pow(r + dr, sigma)));
					integrand /= q * sqrt((q / pow(r + dr, sigma)) * (q / pow(r + dr, sigma)) - 1);

					integrand_next = cosh(l / sigma * acosh((q + dq) / pow(r + dr, sigma)));
					integrand_next /= (q + dq)  * sqrt(((q + dq) / pow(r + dr, sigma)) * ((q + dq) / pow(r + dr, sigma)) - 1);

					qIntegralSumR += 0.5 * (integrand * fhatR + integrand_next * fhatR_next) * dq;
					qIntegralSumG += 0.5 * (integrand * fhatG + integrand_next * fhatG_next) * dq;
					qIntegralSumB += 0.5 * (integrand * fhatB + integrand_next * fhatB_next) * dq;

					q += dq;
				}

				// approx. derivative of the q-integral \w respect to r
				qIntegralR_dr = (qIntegralSumR - qIntegralR_dr) / dr;
				qIntegralG_dr = (qIntegralSumG - qIntegralG_dr) / dr;
				qIntegralB_dr = (qIntegralSumB - qIntegralB_dr) / dr;

				// add circular fourier l-term
				cHarmonicSeriesSumR += cos(l * theta) * qIntegralR_dr;
				cHarmonicSeriesSumG += cos(l * theta) * qIntegralG_dr;
				cHarmonicSeriesSumB += cos(l * theta) * qIntegralB_dr;
			}

			cHarmonicSeriesSumR *= seriesRconst;
			cHarmonicSeriesSumG *= seriesRconst;
			cHarmonicSeriesSumB *= seriesRconst;

			if (cHarmonicSeriesSumR < minR) {
				minR = cHarmonicSeriesSumR;
			}
			if (cHarmonicSeriesSumR > maxR) {
				maxR = cHarmonicSeriesSumR;
			}

			if (cHarmonicSeriesSumG < minG) {
				minG = cHarmonicSeriesSumG;
			}
			if (cHarmonicSeriesSumG > maxG) {
				maxG = cHarmonicSeriesSumG;
			}

			if (cHarmonicSeriesSumB < minB) {
				minB = cHarmonicSeriesSumB;
			}
			if (cHarmonicSeriesSumB > maxB) {
				maxB = cHarmonicSeriesSumB;
			}

			resValuesR[i * h + j] = cHarmonicSeriesSumR;
			resValuesG[i * h + j] = cHarmonicSeriesSumG;
			resValuesB[i * h + j] = cHarmonicSeriesSumB;

			// std::cout << "(" << r << ", " << theta << ") --> " << cHarmonicSeriesSumR << std::endl;

			theta += dTheta;
			// end theta cycle
		}
		p += dp;
		r += dr;
		// end r cycle

		if (j % percent == 0) {
			std::cout << "\rInverse Cormack Transform: "
				<< (int)(((double)j / (double)w * 100) + 0.5) << " % complete ";
		}
	}

	std::cout << std::endl;
	std::cout << "RGB limits:" << std::endl;
	std::cout << "R: [" << minR << ", " << maxR << "]" << std::endl;
	std::cout << "G: [" << minG << ", " << maxG << "]" << std::endl;
	std::cout << "B: [" << minB << ", " << maxB << "]" << std::endl;

	// writing values into the result image
	Magick::Image result(Magick::Geometry(w, h), "black");
	result.type(Magick::TrueColorType);
	result.modifyImage();

	Magick::Pixels view(result);
	Magick::Quantum *pixels = view.get(0, 0, w, h);
	Magick::ColorRGB col;

	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			col = Magick::ColorRGB(
				resValuesR[i * h + j] / (maxR - minR) - minR,
				resValuesG[i * h + j] / (maxG - minG) - minG,
				resValuesB[i * h + j] / (maxB - minB) - minB
			);

			*pixels++ = col.quantumRed();
			*pixels++ = col.quantumGreen();
			*pixels++ = col.quantumBlue();
		}

		if (i % percent == 0) {
			std::cout << "\rWriting results : " << (int)(((double)i / (double)h * 100) + 0.5) << " % complete ";
		}
	}
	std::cout << std::endl;
	std::cout << "Inverse Cormack " << paramString << "-transform complete." << std::endl;

	view.sync();

	delete[] resValuesR;
	delete[] resValuesG;
	delete[] resValuesB;

	return result;
}

int main(int /*argc*/, char **argv)
{
	// Initialize ImageMagick install location for Windows
	Magick::InitializeMagick(*argv);

	Magick::Image orig_image;

	try {
		orig_image.read("Shepp-Logan-phantom.pgm");

		double s = 0.3;
		orig_image.resize(Magick::Geometry(s * orig_image.columns(), s * orig_image.rows()));

		// alpha transforms
		Magick::Image result_1 = CormackTransform(&orig_image, 1.);
		//Magick::Image result_2 = CormackTransform(&orig_image, 0.5);

		// beta transforms
		//Magick::Image result_3 = CormackTransform(&orig_image, 1., false);
		//Magick::Image result_4 = CormackTransform(&orig_image, 0.5, false);

		result_1.write("Sinogram1.jpg");
		//result_2.write("Sinogram2.jpg");
		//result_3.write("Sinogram3.jpg");
		//result_4.write("Sinogram4.jpg");

		Magick::Image result_5 = InverseCormackTransform(&result_1, 1.);

		result_5.write("FFT_Sinogram1.jpg");
	}
	catch (Magick::Exception &error_)
	{
		std::cout << "Caught exception: " << error_.what() << std::endl;
		return 1;
	}
	return 0;
}