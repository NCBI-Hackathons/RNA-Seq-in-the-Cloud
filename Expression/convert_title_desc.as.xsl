<xsl:stylesheet version='1.0'
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    <xsl:output omit-xml-declaration="yes" method='text'/>
    <xsl:variable name="tab"><xsl:text>&#x09;</xsl:text></xsl:variable>
    <xsl:variable name="newline"><xsl:text>&#xa;</xsl:text></xsl:variable>


	<xsl:template match="/">
			<xsl:for-each select="/EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/RUN_SET/RUN">
        		<xsl:value-of select="../../EXPERIMENT/STUDY_REF/@accession"/>
        		<xsl:value-of select="$tab" />
			<!-- fetch GEO name if it exists, other source refnames occasionally present too -->
        		<!--xsl:value-of select="../../EXPERIMENT/STUDY_REF/@refname"/-->
        		<!--xsl:value-of select="$tab" /-->
        		<xsl:value-of select="../../EXPERIMENT/@accession"/>
        		<xsl:value-of select="$tab" />
       			<!-- renamed variable xsl:value-of select="@acc"/ -->
       			<xsl:value-of select="@accession"/>
        		<xsl:value-of select="$tab" />
       			<xsl:value-of select="@total_bases"/>
        		<xsl:value-of select="$tab" />
       			<xsl:value-of select="@total_spots"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/@accession"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_LINKS/SAMPLE_LINK/XREF_LINK[DB='biosample']/LABEL"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_NAME/TAXON_ID"/>
        		<xsl:value-of select="$tab" />
        		<xsl:if test="../../EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED">
        		<xsl:text>PAIRED</xsl:text>
        		</xsl:if>
        		<xsl:value-of select="$tab" />
			<xsl:if test="../../EXPERIMENT/PLATFORM">
			<xsl:value-of select="name(../../EXPERIMENT/PLATFORM/*)" />	
        		</xsl:if>
        		<xsl:value-of select="$tab" />
        		<xsl:if test="../../EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE">
        		<xsl:value-of select="../../EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE" />
        		</xsl:if>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='tissue']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='tissue source']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='OrganismPart']/VALUE"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='cell line']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='source cell line']/VALUE"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='developmental stage']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='DevelopmentalStage']/VALUE"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_NAME/ANONYMIZED_NAME"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_NAME/COMMON_NAME"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="substring(../../SAMPLE/TITLE,1,40)"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_CONSTRUCTION_PROTOCOL" />
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../STUDY/DESCRIPTOR/STUDY_TITLE" />
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='strain']/VALUE"/>
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='developmental stage']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='DevelopmentalStage']/VALUE"/>
        		<xsl:value-of select="../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='diseasestatus']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='disease']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='study disease']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='disease state']/VALUE  " />
        		<xsl:value-of select="$tab" />
        		<xsl:value-of select="../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='celltype']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='cell type']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='cell_type']/VALUE | ../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='cell_subtype']/VALUE  " />
        		<xsl:value-of select="$tab" />
                        <xsl:choose>
                        <xsl:when  test="../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='BioSourceProvider']/VALUE">
                        <xsl:value-of  select="../../SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='BioSourceProvider']/VALUE"/>
                        </xsl:when>
                        <xsl:otherwise>
        		<xsl:value-of select="../../STUDY/DESCRIPTOR/STUDY_DESCRIPTION" />
                        </xsl:otherwise>
                        </xsl:choose>
        		<xsl:value-of select="$newline" />
			</xsl:for-each>
	</xsl:template>

</xsl:stylesheet>

