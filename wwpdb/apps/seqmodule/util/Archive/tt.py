
import re, sys, traceback

def getStrain(sourceName=''):
    myRegex=r'\((.*)\)'
    try:
        #self.__lfh.write("+SequenceFeature.decodeUniProtSourceName() decoding %r\n" % sourceName)
        m =re.findall(myRegex, sourceName)
        if len(m)> 0:
            newSourceName = re.sub(myRegex, '', sourceName)
            if str(m[0]).startswith("strain"):
                newStrain=str(m[0])[6:]
            else:
                newStrain=str(m[0])

            if ((len(newStrain) < 1) or (newStrain.upper() == 'NONE')):
                newStrain == ''
            return str(newSourceName).strip(),str(newStrain).strip()
    except:
        sys.stderr.write("+SequenceFeature.decodeUniProtSourceName() failed \n")
        traceback.print_exc(file=self.__lfh)

    return sourceName,''

def xgetStrain(sourceName=''):
    myRegex=r'\((.*)\)'
    try:
        m =re.findall(myRegex, sourceName)
        if len(m)> 0:
            newSourceName = str(re.sub(myRegex, '', sourceName)).strip()
            if str(m[0]).startswith("strain"):
                newStrain=str(m[0])[6:]
            else:
                newStrain=str(m[0])                
            return newSourceName,newStrain.strip()
    except:
        sys.stderr.write("+SequenceFeature.decodeUniProtSourceName() failed \n")
        traceback.print_exc(file=sys.stderr)

    return sourceName,''

if __name__ == '__main__':

    a='Clostridium cellulolyticum (strain ATCC 35319 / DSM 5812 / JCM 6584 / H10)'
    print a 
    s1,s2=getStrain(sourceName=a)
    print s1
    print s2
    #
    a='Clostridium cellulolyticum (strain ATCC 35319 / DSM (5812) / JCM (6584) / (H10)())()'
    print a
    s1,s2=getStrain(sourceName=a)
    print s1
    print s2
