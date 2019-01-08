class Solution:
    def findMedianSortedArrays(self, nums1, nums2):
        """
        :type nums1: List[int]
        :type nums2: List[int]
        :rtype: float
        """
        num = nums1+nums2
        num.sort()
        size = len(num)
        if size%2 == 0 :
            div = int(size/2)
            return (num[div -1]+num[div])/2
        else:
            div = int(size/2)
            return num[div]