#ifndef __itkMultiLocalCreasenessImageFilter_txx
#define __itkMultiLocalCreasenessImageFilter_txx

#include "itkMultiLocalCreasenessImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkStructureTensorImageFilter.h"
#include "itkSymmetricEigenVectorAnalysisImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{

/**
 * Constructor
 */
template <typename TInputImage, typename TOutputImage >
MultiLocalCreasenessImageFilter<TInputImage,TOutputImage>
::MultiLocalCreasenessImageFilter()
{
  // emptyness constructroness
  m_Sigma = 1.0;
}

/**
 * Set value of Sigma
 */
template <typename TInputImage, typename TOutputImage >
void
MultiLocalCreasenessImageFilter<TInputImage,TOutputImage>
::SetSigma( RealType sigma )
{

  m_Sigma = sigma;

  this->Modified();

}

/**
 * Set value of Rho
 */
template <typename TInputImage, typename TOutputImage >
void
MultiLocalCreasenessImageFilter<TInputImage,TOutputImage>
::SetRho( RealType rho )
{

  m_Rho = rho;

  this->Modified();

}

template <typename TInputImage, typename TOutputImage>
void
MultiLocalCreasenessImageFilter<TInputImage,TOutputImage>
::GenerateInputRequestedRegion() throw(InvalidRequestedRegionError)
{
  // call the superclass' implementation of this method. this should
  // copy the output requested region to the input requested region
  Superclass::GenerateInputRequestedRegion();

  // This filter needs all of the input
  typename MultiLocalCreasenessImageFilter<TInputImage,TOutputImage>
    ::InputImagePointer image = 
      const_cast<InputImageType *>( this->GetInput() );
  image->SetRequestedRegion( this->GetInput()->GetLargestPossibleRegion() );
}

template <typename TInputImage, typename TOutputImage>
void
MultiLocalCreasenessImageFilter<TInputImage,TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *output)
{
  TOutputImage *out = dynamic_cast<TOutputImage*>(output);

  if (out)
    {
    out->SetRequestedRegion( out->GetLargestPossibleRegion() );
    }
}


template <typename TInputImage, typename TOutputImage >
void
MultiLocalCreasenessImageFilter<TInputImage,TOutputImage >
::GenerateData(void)
{

  // Create a process accumulator for tracking the progress of this
  // minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  const typename TInputImage::ConstPointer   inputImage( this->GetInput() );

  typedef itk::StructureTensorImageFilter < InputImageType > 
    StructureTensorFilterType;
  typedef StructureTensorFilterType::OutputImageType 
    SymmetricSecondRankTensorImageType;
  StructureTensorFilterType::Pointer strucTensFilter = 
    StructureTensorFilterType::New();
  strucTensFilter->SetInput( inputImage );
  strucTensFilter->SetSigma(m_Sigma);
  strucTensFilter->SetRho(m_Rho);
  //strucTensFilter->ReleaseDataFlagOn();

  typedef  itk::FixedArray< InternalRealType, ImageDimension> 
    EigenValueArrayType;
  typedef  itk::Image< EigenValueArrayType, ImageDimension> EigenValueImageType;

  // Generate eigen vector image
  typedef  itk::Matrix< InternalRealType, ImageDimension, ImageDimension> 
    EigenVectorMatrixType;
  typedef  itk::Image< EigenVectorMatrixType, ImageDimension> 
    EigenVectorImageType;
  typedef itk::
    SymmetricEigenVectorAnalysisImageFilter<SymmetricSecondRankTensorImageType,
      EigenValueImageType, EigenVectorImageType> EigenVectorAnalysisFilterType;
  EigenVectorAnalysisFilterType::Pointer eigenVectorAnalysisFilter = 
    EigenVectorAnalysisFilterType::New();
  eigenVectorAnalysisFilter->SetDimension( ImageDimension );
  eigenVectorAnalysisFilter->OrderEigenValuesBy( 
		EigenVectorAnalysisFilterType::FunctorType::OrderByMagnitude);
  eigenVectorAnalysisFilter->SetInput( strucTensFilter->GetOutput() );   
  
  //Generate an image with eigen vector pixel that correspond to the largest
  // eigen value 
  EigenVectorImageType::ConstPointer eigenVectorImage =
    eigenVectorAnalysisFilter->GetOutput();
  
  typedef itk::Image< CovariantVector< InternalRealType, ImageDimension >,
    ImageDimension > GradientImageType;
  
  // Compute the contribution of each filter to the total progress.
  const double weight = 1.0 / 2.0;
  progress->RegisterInternalFilter( strucTensFilter, weight );
  progress->RegisterInternalFilter( eigenVectorAnalysisFilter, weight);
  progress->ResetProgress();

  eigenVectorAnalysisFilter->Update();

  // reorient eigenvectors according to gradient.  
  // first, create a image to store the primary eigenvectors
  typedef itk::Image< CovariantVector< InternalRealType, ImageDimension >,
    ImageDimension > VectorImageType;
  VectorImageType::Pointer primaryEigenVectorImage = VectorImageType::New();

  VectorImageType::RegionType region;
  region.SetSize(eigenVectorImage->GetLargestPossibleRegion().GetSize());
  region.SetIndex(eigenVectorImage->GetLargestPossibleRegion().GetIndex());
  primaryEigenVectorImage->SetRegions( region );
  primaryEigenVectorImage->SetOrigin(eigenVectorImage->GetOrigin());
  primaryEigenVectorImage->SetSpacing(eigenVectorImage->GetSpacing());
  primaryEigenVectorImage->Allocate();

  //Fill up the buffer with null vector
  itk::CovariantVector< InternalRealType, ImageDimension > nullVector
    ( ImageDimension );
  for ( unsigned int i=0; i < ImageDimension; i++ )
    {
    nullVector[i] = 0.0;
    }
  primaryEigenVectorImage->FillBuffer( nullVector );
  
  itk::ImageRegionConstIterator<GradientImageType> gradIt(
    strucTensFilter->GetSmoothGradient(),
    strucTensFilter->GetSmoothGradient()->GetLargestPossibleRegion());
  itk::ImageRegionIterator<EigenVectorImageType> eigenVecIt(
    eigenVectorAnalysisFilter->GetOutput(),
    eigenVectorAnalysisFilter->GetOutput()->GetLargestPossibleRegion());
  itk::ImageRegionIterator<VectorImageType> outeigenVectIt(
    primaryEigenVectorImage,
    primaryEigenVectorImage->GetLargestPossibleRegion());
  
  // Important: each ROW of the Matrix (EigenVectorImageType) is a eigenvector.
  // and the last row is the most significant eigenvector.
  for ( gradIt.GoToBegin(), eigenVecIt.GoToBegin(), outeigenVectIt.GoToBegin();
    !gradIt.IsAtEnd(); ++gradIt, ++eigenVecIt, ++outeigenVectIt) 
    {
    GradientImageType::PixelType grad = gradIt.Get();
    EigenVectorImageType::PixelType eigen = eigenVecIt.Get();
    double sign = 0;
    for (size_t idx = 0; idx < ImageDimension; ++idx) 
      {
      sign += grad[idx]*eigen[ImageDimension-1][idx];
      }
    if(sign < 0) 
      {
      for (size_t idx = 0; idx < ImageDimension; ++idx)
        {
        eigen[ImageDimension-1][idx] = - eigen[ImageDimension-1][idx];
        }
      } 
    else if( sign == 0) 
      {
      for (size_t idx = 0; idx < ImageDimension; ++idx) 
        {
        eigen[ImageDimension-1][idx] = 0;
        }
      }
    VectorImageType::PixelType vect;
    for (size_t idx = 0; idx < ImageDimension; ++idx) 
      {
      vect[idx] = eigen[ImageDimension-1][idx];
      }
    outeigenVectIt.Set(vect);
    }
  // primaryEigenVectorImage now contains the primary eigenvector correctly 
  // oriented. Now compute gradients of reoriented eigenvalues;
  m_outputImage = this->GetOutput();
  m_outputImage->SetLargestPossibleRegion(
    inputImage->GetLargestPossibleRegion() );
  m_outputImage->SetBufferedRegion(
    inputImage->GetBufferedRegion() );
  m_outputImage->SetRequestedRegion(
    inputImage->GetRequestedRegion() );
  m_outputImage->Allocate();
  m_outputImage->FillBuffer(0);
  
  itk::ImageRegionIterator<OutputImageType> oit(
    m_outputImage, m_outputImage->GetLargestPossibleRegion());
    
  typedef NthElementImageAdaptor< VectorImageType, InternalRealType > 
    VectorImageAdaptorType;
  VectorImageAdaptorType::Pointer nthAdaptor = VectorImageAdaptorType::New();
  nthAdaptor->SetImage(primaryEigenVectorImage);
  
  typedef itk::GradientRecursiveGaussianImageFilter<VectorImageAdaptorType,
    GradientImageType> GradientAdaptorFilter;
    GradientAdaptorFilter::Pointer gradientAdaptorFilter = 
      GradientAdaptorFilter::New();
  for (size_t idx = 0; idx < ImageDimension; ++idx) 
    {
    nthAdaptor->SelectNthElement( idx );    
    gradientAdaptorFilter->SetInput(nthAdaptor);
    gradientAdaptorFilter->SetSigma(m_Sigma);
    gradientAdaptorFilter->Update();    
    gradientAdaptorFilter->ReleaseDataFlagOn();
    itk::ImageRegionConstIterator<GradientImageType> git(
      gradientAdaptorFilter->GetOutput(),
      gradientAdaptorFilter->GetOutput()->GetLargestPossibleRegion());
    for ( oit.GoToBegin(), git.GoToBegin(); !oit.IsAtEnd(); ++oit, ++git ) 
      {
      GradientImageType::PixelType data = git.Get();
      InternalRealType tmp = oit.Get();
      oit.Set(tmp - data[idx]);    
      }
    }
}

template <typename TInputImage, typename TOutputImage>
void
MultiLocalCreasenessImageFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << "NormalizeAcrossScale: " << m_NormalizeAcrossScale << std::endl;
}

} // end namespace itk

#endif
