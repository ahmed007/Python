def intersection(self, nums1, nums2):
        """
        :type nums1: List[int]
        :type nums2: List[int]
        :rtype: List[int]
        """
        
       a=set()
       b=set()
       c=set()
       for i in nums1:
            a.add(i)
       for v in nums2:
            b.add(v)
    
       c=a.intersection(b)
        
       return list(c)
	   