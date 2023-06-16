module helpers
   implicit none

contains
   ! the `numbers` array MUST be sorted. It returns the index of `searchVal` in numbers
   ! if it exists or -1 if it doesn't
   function binarySearch(numbers, searchVal) result(index)
      integer, allocatable, intent(in) :: numbers(:)
      integer, intent(in) :: searchVal
      integer :: index, minIndex, midIndex, maxIndex

      if(.not. allocated(numbers) .or. size(numbers) == 0) return
      index = -1
      ! fortran does integer division automatically if the two operands are integers
      minIndex = 1
      maxIndex = size(numbers)
      do while (minIndex <= maxIndex)
         midIndex = minIndex + (maxIndex - minIndex)/2
         if ( numbers(midIndex) == searchVal ) then
            index = midIndex
            exit
         else if ( numbers(midIndex) < searchVal ) then
            ! shift min
            minIndex = midIndex+1
         else
            ! shift max: numbers(midIndex)>searchVal
            maxIndex = midIndex-1
         end if
      end do
   end function binarySearch
end module helpers
