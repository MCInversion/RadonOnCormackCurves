#include "pch.h"
#include <iostream>
#include <Magick++.h>
#include <cmath>
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

Magick::Image CormackTransformAlpha(Magick::Image *image, double alpha) {
	const double phiMin = 0.; const double phiMax = 2 * M_PI;
	const double pMin = 0.; const double pMax = 1.;
	const double psiMax = 0.99 * M_PI / (2 * alpha);
	const double psiMin = -psiMax;
	double phi, p, dArcLength;
	double fvalR, fvalG, fvalB; // , fvalR_prev, fvalG_prev, fvalB_prev;
	double psi, r, theta, x, y, psiIntegSumR, psiIntegSumG, psiIntegSumB;
	int imgX, imgY;
	int imgP, imgPhi;
	Magick::ColorRGB px;

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

	std::cout << "initializing Cormack alpha transform \\w alpha = " << alpha << std::endl;
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

				r = p / pow(cos(alpha * psi), 1. / alpha);
				theta = phi + psi;

				x = polarToX(r, theta);
				y = polarToY(r, theta);

				imgX = XtoPixelX(w, x);
				imgY = YtoPixelY(h, y);

				if (!isInsideImage(imgX, imgY, h, w)) {
					continue;
				}

				px = image->pixelColor(imgX, imgY);

				fvalR = px.red(); fvalG = px.green(); fvalB = px.blue();

				// fvalR = fvalG = fvalB = f_gaussians(x, y);

				dArcLength = dPsi / pow(cos(alpha * psi), (1. + 1. / alpha));

				psiIntegSumR += p * fvalR * dArcLength;
				psiIntegSumG += p * fvalG * dArcLength;
				psiIntegSumB += p * fvalB * dArcLength;
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
	std::cout << "Cormack alpha-transform complete." << std::endl;


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

		double s = 0.5;
		orig_image.resize(Magick::Geometry(s * orig_image.columns(), s * orig_image.rows()));

		Magick::Image result_1 = CormackTransformAlpha(&orig_image, 1.);
		Magick::Image result_2 = CormackTransformAlpha(&orig_image, 0.5);

		result_1.write("Sinogram1.jpg");
		result_2.write("Sinogram2.jpg");
	}
	catch (Magick::Exception &error_)
	{
		std::cout << "Caught exception: " << error_.what() << std::endl;
		return 1;
	}
	return 0;
}