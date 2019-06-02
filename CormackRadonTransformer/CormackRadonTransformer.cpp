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

	std::complex<double> *outR = new std::complex<double>[w];
	std::complex<double> *outG = new std::complex<double>[w];
	std::complex<double> *outB = new std::complex<double>[w];

	double minR = DBL_MAX;
	double maxR = -DBL_MAX;
	double minG = minR;
	double maxG = maxR;
	double minB = minR;
	double maxB = maxR;

	for (int i = 0; i < h; i++) {
		// fill 3 1D vectors
		for (int j = 0; j < w; j++) {
			px = image->pixelColor(i, j);

			rowR[j] = px.red();
			rowG[j] = px.green();
			rowB[j] = px.blue();
		}
		fft(rowR, outR, w);
		fft(rowG, outG, w);
		fft(rowB, outB, w);

		for (int j = 0; j < w; j++) {
			if (outR[j].real() < minR) {
				minR = outR[j].real();
			}
			if (outR[j].real() > maxR) {
				maxR = outR[j].real();
			}

			if (outG[j].real() < minG) {
				minG = outG[j].real();
			}
			if (outG[j].real() > maxG) {
				maxG = outG[j].real();
			}

			if (outB[j].real() < minB) {
				minB = outB[j].real();
			}
			if (outB[j].real() > maxB) {
				maxB = outB[j].real();
			}
		}

		for (int j = 0; j < w; j++) {
			col = Magick::ColorRGB(
				outR[j].real() / (maxR - minR) - minR,
				outG[j].real() / (maxG - minG) - minG,
				outB[j].real() / (maxB - minB) - minB
			);

			*pixels++ = col.quantumRed();
			*pixels++ = col.quantumGreen();
			*pixels++ = col.quantumBlue();
		}

		if (i % percent == 0) {
			std::cout << "\rFFT : " << (int)(((double)i / (double)h * 100) + 0.5) << " % complete ";
		}
	}
	view.sync();
	std::cout << std::endl;

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
	double dTheta = 2 * M_PI / integralSize;
	int percent = (h > 100 ? (int)round((double)h / 100.) : 1);
	int N = 5;

	double qIntegralSumR, qIntegralSumG, qIntegralSumB;
	double fhatR, fhatG, fhatB, fhatR_next, fhatG_next, fhatB_next;
	double qMin, qMax = 2. * std::max(h, w);

	std::string paramString = alpha_transform ? "alpha" : "beta";
	std::cout << "initializing Inverse Cormack " << paramString << "-transform \\w " << paramString << " = " << sigma << ", Fourier order: " << N << std::endl;
	std::cout << "image : " << h << "x" << w << std::endl;

	for (int j = 0; j < w; j++) {

		for (int i = 0; i < h; i++) {

			qIntegralSumR = 0.; qIntegralSumG = 0.; qIntegralSumB = 0.;
			for (int q = 0; q < integralSize; q++) {


			}

			if (qIntegralSumR < minR) {
				minR = qIntegralSumR;
			}
			if (qIntegralSumR > maxR) {
				maxR = qIntegralSumR;
			}

			if (qIntegralSumG < minG) {
				minG = qIntegralSumG;
			}
			if (qIntegralSumG > maxG) {
				maxG = qIntegralSumG;
			}

			if (qIntegralSumB < minB) {
				minB = qIntegralSumB;
			}
			if (qIntegralSumB > maxB) {
				maxB = qIntegralSumB;
			}
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

	return fourier_img;
}

int main(int /*argc*/, char **argv)
{
	// Initialize ImageMagick install location for Windows
	Magick::InitializeMagick(*argv);

	Magick::Image orig_image;

	try {
		orig_image.read("Shepp-Logan-phantom.pgm");

		double s = 0.5;
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

		result_5.write("IFFT_Sinogram1.jpg");
	}
	catch (Magick::Exception &error_)
	{
		std::cout << "Caught exception: " << error_.what() << std::endl;
		return 1;
	}
	return 0;
}