import json
import sys
import md5
import base64
import re
import binascii
import urllib2


API_URL = 'http://seqdepot.net/api/v1'
VERSION = '0.01'



"""-------------------------------------------------------------------------------------------------------------------
NAME:
SeqDepot - SeqDepot database utility class

SYNOPSIS:
import SeqDepot_p2

DESCRIPTION:
The SeqDepot class facilitates interacting with the public SeqDepot
database via its RESTful interface. This package also provides helper
methods such as a FASTA parser routine for working with sequences,
and methods for manipulating identifiers.

SUBROUTINES/METHODS:
new

Creates a new SeqDepot instance.

"""


class new(object):
    def __init__(self):
        # Internal string buffer used for reading fasta sequences
        self.fastaBuffer = None
        # Store any error here
        self.lastError = None
        # After querying the REST api for tool info, cache the result
        self.toolCache = None

        self.toolPosition = {}


    def aseqIdFromMD5Hex(self, MD5hex):
        """aseqIdFromMD5Hex

Static method that converts an MD5 hexadecimal string into its aseqId
equivalent.

Parameters:
  MD5hex: hexadecimal MD5 string

Return:
  aseqId

"""

        return base64.encodestring(binascii.unhexlify(MD5hex)).replace('/','_').replace('=','').replace('+','-').replace('\n','')

    def aseqIdFromSequence(self, sequence):
        """aseqIdFromSequence

Static method for computing the aseqId for a given sequence. It is
recommended that all sequences are cleaned before calling this method.

Parameters:
  sequence

Returns:
  aseqId
"""
        return base64.encodestring(md5.new(sequence.replace('-','')).digest()).replace('/','_').replace('=','').replace('+','-').replace('\n','')

    def cleanSequence(self, sequence):
        """cleanSequence

Static method for removing all whitespace characters from sequence
and replaces all digits or non-word characters with an ampersand
character (for easy identification of invalid symbols).

Parameters:
  sequence: sequence character string

Returns:
  string
"""

        if sequence.__class__ != str:
            print 'Missing sequence parameter'
        else:
            sequence = re.sub(r'\s','', sequence)
            sequence = re.sub(r'\W|\d','@',sequence)
        return sequence

    def find(self, ids, params={}):
        """ find

Retrieves fields for aseqIds from SeqDepot. All fields are returned
by default; however, a subset may be returned by specifying the desired
fields (see below).

Returns a mixed array of hashes or nulls, indicating whether the
respective requested aseqId was found. Any nulls in the array indicate
that the requested Aseq ID was not found - not that some other error
occurred.

Parameters:
  ids: single or array of Aseq ID | GI | UniProt ID | PDB ID | MD5 hex

  dictionary: accepts the following fields

      type: <string> type of identifier defaults to aseq_id, but may also be
        gi, uni, pdb, md5_hex

      fields (optional): comma-separated string of primary fields, secondary
        fields must be listed in parentheses, separated with pipe (|)
        symbols, and immediately suffix the respective primary field name.
        For example:

          s,l   --> Returns the sequence and length primary fields

          l,t,x --> Returns the length and all tool and cross-references

          l,t(pfam26|smart) --> Returns the length primary field and pfam and
            smart secondary fields of t

          l,t(pfam26|smart),x(gi) --> Same as the above plus any gi numbers

      labelToolData: <boolean> defaults to False; if True, converts any tool
        data (the t field) into an array of dictionary with meaningful names.

Returns:
  (null | mixed array of dictionary and nulls): dictionary if the given Aseq ID was
    successfully located or null otherwise; null indicates that an
    error occurred - call lastError() for details.
"""
        if ids.__class__ == list:
            ids_list = [str(i) for i in ids]
        elif ids.__class__ == str or ids.__class__ == int:
            ids_list = [ids]
        else:
            print 'Missing ids parameter or invalid id type (string or list)'

        self.clearError_()

        url = API_URL + '/aseqs'
        stringyParams = self.urlParams_(params)

        url += '?' + stringyParams
        request = urllib2.Request(url,'\n' + '\n'.join(ids_list))
        content = self.lwpResponse(request)

        results = []
        re_iter = re.finditer(r'(\S.*)', content.read())
        for i in re_iter:
            line = i.groups()[0]
            queryId, code, aseID, json_f = line.split('\t')[:4]
            try:
                code = int(code)
            except ValueError:
                code = None
            result = { 'code':code, 'query':queryId}
            if code == 200:
                result['data'] = json.loads(json_f)
                if 'labelTooldata' in params.keys() and 't' in result['data'].keys():
                    result['data']['t'] = self.labelToolData_(result['data']['t'])
            results.append(result)
        return results

    def findOne(self, ids, params={}):
        """findOne

Retrieves fields for aseqId from SeqDepot. All fields are returned
by default; however, a subset may be returned by specifying the desired
fields (see find for details).

A null return value indicates that either the aseqId does not exist
in SeqDepot or an error occurred. Call lastError() to retrieve
details about any errors. If the aseqId does not exist, then
lastError() will return null.

Parameters:
  Aseq ID | GI | UniProt ID | PDB ID | MD5 hex
  params (optional, see find documentation)

Returns:
  Dictionary if successfully found; null if not found or an error occurred -
     call lastError() for details.
"""
        self.clearError_()
        ids = str(ids)
        url = API_URL + '/aseqs/' + ids
        stringyParams = self.urlParams_(params)
        url += '?' + stringyParams
        content = self.lwpResponse(url)
        result = json.loads(content.read())
        if 'labelTooldata' in params.keys() and 't' in result.keys():
            result['t'] = self.labelToolData_(result['t'])
        return result

    def isToolDone(self, toolId = None, status = None):
        """isToolDone

Returns True if the requested tool has been marked as done from the
status string. The status string corresponds to the _s field in
SeqDepot and contains information about which predictive self.toolCache have
been executed and whether any results were found with this tool. See
the main documentation for more details.

Parameters:
  toolId <string>
  status <string>

Returns:
  boolean

Null is returned if toolId is not valid.
"""
        if toolId == None:
            return False
        if status == None:
            return False

        tools = self.tools()
        if not self.toolCache:
            print "Unable to fetch tool metadata from SeqDepot: " + self.lastError + "\n";

#############PATCH - Handling wrong ToolID ##########################
        if toolId in self.toolPosition.keys():
            pos = self.toolPosition[toolId]
        else:
            print "ToolId not recognized"
            return False
#################END PATCH ##########################################
#################ORIGINAL############################################
#        pos = self.toolPosition[toolId]
#####################################################################
        if pos.__class__ is int:
            if pos < len(status):
                if status[pos] != '-':
                    return True
        return False

    def isValidAseqId(self, aseqId = None):
        """isValidAseqId

Static method that returns True if aseqId is validly formatted; False otherwise
"""
        if aseqId == None:
            return False
        if re.match('[A-Za-z0-9_-]{22}$',aseqId):
            return True
        else:
            return False

    def isValidFieldString(self, fields=None):
        """ isValidFieldString

Checks if the requested field string is valid.

Parameters:
  fields

Returns:
  boolean
"""
        if not fields:
            return True

        primaries = str(fields).split(',')
        for primary in primaries:
            p = re.match('^([A-Za-z_-][A-Za-z0-9-_]*)(?:\(([A-Za-z0-9-_|]+)\))?$', primary)
            if not p:
                return False

            if p.groups()[1]:
                for subField in p.groups()[1].split('|'):
                    if not subField:
                        return False
        return True

    def lastError(self):
        """lastError

Returns any error that may have occurred or null if there was no error
for the last find operation.
"""
        return self.lastError


    def MD5HexFromAseqId(self, aseqId = ''):
        """ MD5HexFromAseqId

Static method that converts an aseqId to its equivalent MD5 hexadecimal
representation.

Parameters:
  aseqId

Returns:
  MD5 hexadecimal string
"""
        if aseqId == '':
            print "Missing aseqId parameter"

        if not self.isValidAseqId(aseqId):
            print "Converting invalid Aseq ID: " + aseqId
        aseqId = aseqId.replace('-', '+').replace('_', '/')
        aseqId += '=' * (24 - len(aseqId))
        return    binascii.hexlify(base64.decodestring(aseqId))

    def MD5HexFromSequence(self, sequence=''):
        """MD5HexFromSequence

Static method for computing the hexadecimal MD5 digest for a given
sequence. It is recommended that all sequences are cleaned before
calling this method.

Parameters:
  sequence

Returns:
  MD5 hexadecimal string
"""
        if sequence == '':
            print "Missing sequence parameter"
            return None
        return binascii.hexlify(md5.new(sequence).digest())

    def primeFastaBuffer(self, buffer):
        """primeFastaBuffer

Sets the internal fastaBuffer to the argument, fastaBuffer. This is useful
when an input stream has already been partially read but not processed as
part of the FASTA parsing. For example, when reading a line from STDIN to
determine if it is FASTA data.

Parameters:
  fastaBuffer: <string>
"""
        self.fastaBuffer = buffer

    def readFastaSequence(self, fh):
        """readFastaSequence

Reads a FASTA-formatted sequence from an open file handle and returns
an array containing the header and the cleaned sequence. The header
will not contain the > symbol. Returns null if there are no more
sequences to be read from the file.

Whitespace is trimmed from both ends of the header line.

Parameters:
  fh: non-null open filehandle

Return:
  Array containing the header and cleaned sequence if a sequence was
    read; None otherwise

"""
        line = self.fastaBuffer if self.fastaBuffer else fh.readline()
        self.fastaBuffer = None
        while line:
            if re.match('^\s*$', line):
                line = fh.readline()
                continue

            if line[0] != '>':
                raise Exception('Invalid FASTA file. Header line must begin with a greater than symbol\nLine: ' + line + '\n\n')

            line = re.sub('\s+$', '', line)
            header = line[1:]
            header = re.sub('^\s+', '', header).decode('string_escape')
            sequence = ''
            line = fh.readline()
            while line:
                if line[0] != '>':
                    sequence += line
                    line = fh.readline()
                    continue

            # We got the next header line. Save it for the next call to this method.
                self.fastaBuffer = line
                break

            sequence = self.cleanSequence(sequence)
            return [header, sequence]

        return None

    def resetFastaBuffer(self):
        """resetFastaBuffer

Clears the internal buffer used to read FASTA sequences. Call this
method before readFastaSequence if all of the following are true:
1) changing filehandles
2) the filehandle has been partially read from
3) the filehandle has not been completely read through to the end

"""
        self.fastaBuffer = None

    def saveImage(self, idss = '', fileName = None, params = {}):
        """saveImage

Saves an image of the corresponding aseq.

Parameters:
  id: <string|number> Aseq ID | GI | UniProt ID | PDB ID | MD5 hex
  fileName: <string> optional; defaults to the id with the appropriate
    image extension
  params: <dictionary>
    fields => <string> see find documentation
    type => type of id (aseq_id | gi | uni | pdb | md5_hex)
    format => png | svg: type of image to save; defaults to png
"""
        ids = str(idss)
        url = API_URL + '/aseqs/' + ids
        format_s = 'png'
        self.clearError_()
        if params != {}:
            if params.keys() == ['format']:
                if params['format'] == 'svg':
                    format_s = 'svg'
        url += '.' + format_s
        stringyParams = self.urlParams_(params)
        url += '?' + stringyParams
        if fileName is None:
            fileName = ids + '.' + format_s
        response = self.lwpResponse(url)
        if response:
            if format_s == 'png':
                output = open(fileName,'wb')
            else:
                output = open(fileName,'w')
            output.write(response.read())
            output.close()
            return 1
        else:
            return None

    def toolFields(self, toolName):
        """toolFields

Returns the field names associated with toolName.

Parameters:
  toolName: string of toolname

Returns:
  array of strings (empty array if toolName is not valid) OR
  None if an error occurred.

"""
        self.toolCache = self.tools()
        if self.toolCache and toolName in self.toolCache.keys():
            if self.toolCache[toolName]['f']:
                return self.toolCache[toolName]['f']
            else:
                return None
        else:
            return None

    def tools(self):
        """tools()

Returns dictionary of predictive self.toolCache available at SeqDepot
"""
        if self.toolCache != None:
            return self.toolCache
        else:
            content = self.lwpResponse(API_URL + '/tools.json').read()

        orderedTools = json.loads(content)['results']

        # Create a tool position lookup
        for i in range(len(orderedTools)):
            idt = orderedTools[i]['id']
            self.toolPosition[idt] = i
        myhash = { i['id']:i for i in orderedTools}
        self.toolCache = myhash
        return self.toolCache

    def toolNames(self):
        """toolNames

Returns an array of tool names used by SeqDepot.
"""
        return self.tools().keys()



#------------------------------------------------------------------------------------------------------
# Private Methods:

    def lwp_(self):
        return

    def clearError_(self):
        self.lastError = None

    def urlParams_(self, params = {}):
        #suggestion: Since there are only few valid type, why not check prior to url submition?
        par_list = []
        if 'fields' in params.keys():
            if self.isValidFieldString(params['fields']):
                par_list.append('fields=' + params['fields'])
            else:
                print 'Invalid fields parameter'

        if 'type' not in params.keys():
            params['type'] = 'aseq_id'

        par_list.append('type='+params['type'])
        return '&'.join(par_list)

    def labelToolData_(self,t={}):
        result = {}
        for toolId in t.keys():
            rows = t[toolId]
            if rows.__class__ != list:
                result[toolId] = rows
                continue

            fieldNames = self.toolFields(toolId)
            hashes = []
            for i in rows:
                myhash = {}
                myhash[fieldNames] = row
                hashes.append(myhash)
            result[toolId] = hashes
        return result

    def lwpResponse(self,request=None):
        try:
            response = urllib2.urlopen(request)
        except urllib2.URLError, e:
            print '\n\nUnable to connect to server; timeout or other internal error'
            error = json.loads(e.read())
            self.lastError = error['message']
            return None
        return response

"""AUTHOR

Davi Ortega, <davi.ortega at gmail.com>
Luke Ulrich, <ulrich.luke+sci at gmail.com>

BUGS

Please report any bugs or feature requests to the AUTHOR.

LICENSE AND COPYRIGHT

Copyright 2013 Luke Ulrich.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://docs.python.org/2/license.html for more information.

"""



